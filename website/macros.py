import urllib
import json

def define_env(env):
    """Define variables, macros, and filters.

    Args:
        env: The environment object used for defining macros and variables.

    Returns:
        None
    """

    @env.macro
    def include_file(filename: str) -> str:
        """Include the contents of a file.

        Args:
            filename (str): The name of the file to include.

        Returns:
            str: The content of the file with specific replacements made,
                 or an error message if the file is not found.

        Raises:
            FileNotFoundError: If the specified file cannot be found.
        """
        try:
            with open(filename, "r", encoding="utf-8") as file:
                content = file.read().replace("website/docs/", "")
                content = content.replace(
                    "[AUTHORS](AUTHORS)", "[AUTHORS](./authors.md)"
                )
                return content
        except FileNotFoundError:
            return f"File not found: {filename}"

    @env.macro
    def get_project_version() -> str:
        url = "https://gitlab.com/api/v4/projects/SFCGAL%2FSFCGAL/releases?per_page=1"
        try:
            with urllib.request.urlopen(url, timeout=5) as response:
                releases = json.loads(response.read())
                if releases:
                    return releases[0]["tag_name"].lstrip("v")
                return "No release found"
        except Exception as e:
            return f"Error fetching version: {e}"