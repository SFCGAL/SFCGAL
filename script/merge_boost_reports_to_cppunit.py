"""
merge_boost_reports.py
----------------------
This python script has been produced by Claude

Reads multiple test reports in Boost XML format and produces a merged report in CppUnit XML format.

Input format (Boost):
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

    # With glob :
    python merge_boost_reports.py boost_report_*.xml -o merged.xml
"""

import argparse
import glob
import sys
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET


# ---------------------------------------------------------------------------
# Intermediate Data Structures
# ---------------------------------------------------------------------------


@dataclass
class TestCaseData:
    """Represents a test case extracted from a Boost report."""

    suite_name: str  # parent <TestSuite> name
    name: str  # <TestCase> name
    isSuccess: bool  # "passed" or "failed" ==> true or false
    assertions_passed: int
    assertions_failed: int


@dataclass
class MergedData:
    """Reading result and merging all input files."""

    test_cases: list[TestCaseData] = field(default_factory=list)

    # Aggregated statistics
    @property
    def total_tests(self) -> int:
        return len(self.test_cases)

    @property
    def failures(self) -> int:
        """Number of failed test cases."""
        return sum(not tc.isSuccess for tc in self.test_cases)

    @property
    def assertions_failed(self) -> int:
        return sum(tc.assertions_failed for tc in self.test_cases)

    @property
    def successful(self) -> list[TestCaseData]:
        return [tc for tc in self.test_cases if tc.isSuccess]

    @property
    def failed(self) -> list[TestCaseData]:
        return [tc for tc in self.test_cases if not tc.isSuccess]


# ---------------------------------------------------------------------------
# Reading Boost reports
# ---------------------------------------------------------------------------


def parse_boost_report(path: Path) -> list[TestCaseData]:
    """
        Parses a Boost XML report and returns the list of Test Case Data found
    Only the direct <Test Case> children of the second level <Test Suite> are
    collected (we ignore the root sequence).
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        sys.exit(f"Parsing error in {path}: {exc}")

    root = tree.getroot()
    if root.tag != "TestResult":
        sys.exit(f"Unexpected root '{root.tag}' in {path} (expected: TestResult)")

    root_suite = root.find("TestSuite")
    if root_suite is None:
        sys.exit(f"No <TestSuite> found in {path}")

    cases: list[TestCaseData] = []

    # Recursive traversal of all child sequences
    def walk(suite: ET.Element) -> None:
        suite_name = suite.get("name", "")
        for child in suite:
            if child.tag == "TestCase":
                cases.append(
                    TestCaseData(
                        suite_name=suite_name,
                        name=child.get("name", ""),
                        isSuccess=child.get("result", "passed") == "passed",
                        assertions_passed=int(child.get("assertions_passed", 0)),
                        assertions_failed=int(child.get("assertions_failed", 0)),
                    )
                )
            elif child.tag == "TestSuite":
                walk(child)

    for child_suite in root_suite:
        if child_suite.tag == "TestSuite":
            walk(child_suite)

    return cases


def merge_reports(input_files: list[Path]) -> MergedData:
    """Merges data from all Boost files into one MergedData."""
    merged = MergedData()
    for path in input_files:
        merged.test_cases.extend(parse_boost_report(path))
    return merged


# ---------------------------------------------------------------------------
# Generating the CppUnit XML report
# ---------------------------------------------------------------------------


def build_cppunit_xml(data: MergedData) -> ET.Element:
    """Constructs the CppUnit XML tree from the merged data."""
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
    # global ID numbering after failed tests
    for idx, tc in enumerate(data.successful, start=len(data.failed) + 1):
        test_el = ET.SubElement(successful_tests_el, "Test", id=str(idx))
        ET.SubElement(test_el, "Name").text = f"{tc.suite_name}::{tc.name}"

    # ── <Statistics> ────────────────────────────────────────────────────────
    stats = ET.SubElement(test_run, "Statistics")
    ET.SubElement(stats, "Tests").text = str(data.total_tests)
    ET.SubElement(stats, "FailuresTotal").text = str(data.failures)
    ET.SubElement(stats, "Errors").text = "0"  # not available in Boost
    ET.SubElement(stats, "Failures").text = str(data.failures)

    return test_run


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Merges multiple Boost XML reports and produces a report au format Cpp Unit XML."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="FILE",
        help="Boost XML files to merge (globs are accepted).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="merged_report.xml",
        metavar="FILE",
        help="CppUnit output file (default: merged report xml).",
    )
    args = parser.parse_args()

    # Resolution of possible globs
    resolved: list[Path] = []
    for pattern in args.inputs:
        matches = [Path(p) for p in glob.glob(pattern)]
        resolved.extend(sorted(matches) if matches else [Path(pattern)])

    if not resolved:
        sys.exit("No input files found.")

    print(f"Merged {len(resolved)} file(s) :")
    for f in resolved:
        print(f"  - {f}")

    # Read + merge
    data = merge_reports(resolved)

    # XML CppUnit building
    cppunit_root = build_cppunit_xml(data)

    # Writing with indentation
    tree = ET.ElementTree(cppunit_root)
    ET.indent(tree, space="  ")  # Python ≥ 3.9

    output_path = Path(args.output)
    with output_path.open("w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n')
        tree.write(fh, encoding="unicode", xml_declaration=False)
        fh.write("\n")

    print(f"\nCppUnit report written in: {output_path}")

    # Summary
    status = "FAILED" if data.failures else "OK"
    print("\n── Summary ───────────────────────────")
    print(f"  Global status : {status}")
    print(f"  All tests     : {data.total_tests}")
    print(f"  Success       : {len(data.successful)}")
    print(f"  Failures      : {data.failures}")
    print(f"  Assertions ✗  : {data.assertions_failed}")
    print("────────────────────────────────────────")


if __name__ == "__main__":
    main()
