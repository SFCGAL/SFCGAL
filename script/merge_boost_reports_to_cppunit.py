"""
merge_boost_reports.py
----------------------
Lit plusieurs rapports de tests au format Boost XML et produit un rapport
fusionné au format CppUnit XML.

Format d'entrée (Boost) :
    <TestResult>
      <TestSuite name="RootSuite" result="..." assertions_passed="..." ...>
        <TestSuite name="SpecificSuite" ...>
          <TestCase name="testFoo" result="passed|failed" .../>
        </TestSuite>
      </TestSuite>
    </TestResult>

Format de sortie (CppUnit) :
    <?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
    <TestRun>
      <FailedTests>
        <FailedTest id="1">
          <Name>SuiteName::testName</Name>
          <FailureType>Assertion</FailureType>
          <Message>N assertions failed</Message>
        </FailedTest>
      </FailedTests>
      <SuccessfulTests>
        <Test id="2">
          <Name>SuiteName::testName</Name>
        </Test>
      </SuccessfulTests>
      <Statistics>
        <Tests>10</Tests>
        <FailuresTotal>1</FailuresTotal>
        <Errors>0</Errors>
        <Failures>1</Failures>
      </Statistics>
    </TestRun>

Usage :
    python merge_boost_reports.py rapport1.xml rapport2.xml [...] -o merged.xml

    # Avec un glob :
    python merge_boost_reports.py boost_report_*.xml -o merged.xml
"""

import argparse
import glob
import sys
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET


# ---------------------------------------------------------------------------
# Structures de données intermédiaires
# ---------------------------------------------------------------------------

@dataclass
class TestCaseData:
    """Représente un cas de test extrait d'un rapport Boost."""
    suite_name: str          # nom de la <TestSuite> parente
    name: str                # nom du <TestCase>
    result: str              # "passed" ou "failed"
    assertions_passed: int
    assertions_failed: int


@dataclass
class MergedData:
    """Résultat de la lecture et fusion de tous les fichiers d'entrée."""
    test_cases: list[TestCaseData] = field(default_factory=list)

    # Statistiques agrégées
    @property
    def total_tests(self) -> int:
        return len(self.test_cases)

    @property
    def failures(self) -> int:
        """Nombre de cas de test en échec."""
        return sum(1 for tc in self.test_cases if tc.result == "failed")

    @property
    def assertions_failed(self) -> int:
        return sum(tc.assertions_failed for tc in self.test_cases)

    @property
    def successful(self) -> list[TestCaseData]:
        return [tc for tc in self.test_cases if tc.result == "passed"]

    @property
    def failed(self) -> list[TestCaseData]:
        return [tc for tc in self.test_cases if tc.result == "failed"]


# ---------------------------------------------------------------------------
# Lecture des rapports Boost
# ---------------------------------------------------------------------------

def parse_boost_report(path: Path) -> list[TestCaseData]:
    """
    Parse un rapport Boost XML et retourne la liste des TestCaseData trouvés.
    Seuls les <TestCase> directs enfants des <TestSuite> de second niveau sont
    collectés (on ignore la suite racine).
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        sys.exit(f"Erreur de parsing dans {path}: {exc}")

    root = tree.getroot()
    if root.tag != "TestResult":
        sys.exit(f"Racine inattendue '{root.tag}' dans {path} (attendu: TestResult)")

    root_suite = root.find("TestSuite")
    if root_suite is None:
        sys.exit(f"Aucune <TestSuite> trouvée dans {path}")

    cases: list[TestCaseData] = []

    # Parcours récursif de toutes les suites enfants
    def walk(suite: ET.Element) -> None:
        suite_name = suite.get("name", "")
        for child in suite:
            if child.tag == "TestCase":
                cases.append(TestCaseData(
                    suite_name=suite_name,
                    name=child.get("name", ""),
                    result=child.get("result", "passed"),
                    assertions_passed=int(child.get("assertions_passed", 0)),
                    assertions_failed=int(child.get("assertions_failed", 0)),
                ))
            elif child.tag == "TestSuite":
                walk(child)

    for child_suite in root_suite:
        if child_suite.tag == "TestSuite":
            walk(child_suite)

    return cases


def merge_reports(input_files: list[Path]) -> MergedData:
    """Fusionne les données de tous les fichiers Boost en un seul MergedData."""
    merged = MergedData()
    for path in input_files:
        merged.test_cases.extend(parse_boost_report(path))
    return merged


# ---------------------------------------------------------------------------
# Génération du rapport CppUnit XML
# ---------------------------------------------------------------------------

def build_cppunit_xml(data: MergedData) -> ET.Element:
    """Construit l'arbre XML CppUnit à partir des données fusionnées."""
    test_run = ET.Element("TestRun")

    # ── <FailedTests> ───────────────────────────────────────────────────────
    failed_tests_el = ET.SubElement(test_run, "FailedTests")
    for idx, tc in enumerate(data.failed, start=1):
        failed_test = ET.SubElement(failed_tests_el, "FailedTest", id=str(idx))
        ET.SubElement(failed_test, "Name").text = f"{tc.suite_name}::{tc.name}"
        ET.SubElement(failed_test, "FailureType").text = "Assertion"
        msg = (
            f"{tc.assertions_failed} assertion(s) failed, "
            f"{tc.assertions_passed} passed"
        )
        ET.SubElement(failed_test, "Message").text = msg

    # ── <SuccessfulTests> ───────────────────────────────────────────────────
    successful_tests_el = ET.SubElement(test_run, "SuccessfulTests")
    # L'id continue la numérotation globale après les tests en échec
    offset = len(data.failed)
    for idx, tc in enumerate(data.successful, start=offset + 1):
        test_el = ET.SubElement(successful_tests_el, "Test", id=str(idx))
        ET.SubElement(test_el, "Name").text = f"{tc.suite_name}::{tc.name}"

    # ── <Statistics> ────────────────────────────────────────────────────────
    stats = ET.SubElement(test_run, "Statistics")
    ET.SubElement(stats, "Tests").text = str(data.total_tests)
    ET.SubElement(stats, "FailuresTotal").text = str(data.failures)
    ET.SubElement(stats, "Errors").text = "0"       # non disponible dans Boost
    ET.SubElement(stats, "Failures").text = str(data.failures)

    return test_run


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Fusionne plusieurs rapports Boost XML et produit un rapport "
            "au format CppUnit XML."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="FILE",
        help="Fichiers Boost XML à fusionner (les globs sont acceptés).",
    )
    parser.add_argument(
        "-o", "--output",
        default="merged_report.xml",
        metavar="FILE",
        help="Fichier de sortie CppUnit (défaut : merged_report.xml).",
    )
    args = parser.parse_args()

    # Résolution des globs éventuels
    resolved: list[Path] = []
    for pattern in args.inputs:
        matches = [Path(p) for p in glob.glob(pattern)]
        resolved.extend(sorted(matches) if matches else [Path(pattern)])

    if not resolved:
        sys.exit("Aucun fichier d'entrée trouvé.")

    print(f"Fusion de {len(resolved)} fichier(s) :")
    for f in resolved:
        print(f"  - {f}")

    # Lecture + fusion
    data = merge_reports(resolved)

    # Construction du XML CppUnit
    cppunit_root = build_cppunit_xml(data)

    # Écriture avec indentation
    tree = ET.ElementTree(cppunit_root)
    ET.indent(tree, space="  ")   # Python ≥ 3.9

    output_path = Path(args.output)
    with output_path.open("w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n')
        tree.write(fh, encoding="unicode", xml_declaration=False)
        fh.write("\n")

    print(f"\nRapport CppUnit écrit dans : {output_path}")

    # Résumé
    status = "FAILED" if data.failures else "OK"
    print("\n── Résumé ──────────────────────────────")
    print(f"  Statut global  : {status}")
    print(f"  Tests total    : {data.total_tests}")
    print(f"  Tests réussis  : {len(data.successful)}")
    print(f"  Tests échoués  : {data.failures}")
    print(f"  Assertions ✗   : {data.assertions_failed}")
    print("────────────────────────────────────────")


if __name__ == "__main__":
    main()
