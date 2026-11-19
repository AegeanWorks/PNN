"""Regression shield for the PNN research documents (PRINCIPLE-005, PRINCIPLE-014)."""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_DOCUMENTS = (
    "README.md",
    "literature_review.md",
    "research_framework.md",
    "agent_design/target_profiles.md",
    "agent_design/bbb_penetration_rules.md",
    "agent_design/oral_bioavailability.md",
    "computational/docking_protocols.md",
    "experimental/experimental_protocols.md",
    "docs/cold_start_guide.md",
    "docs/reproducibility.md",
)


class ProjectDocsTest(unittest.TestCase):
    def test_required_documents_exist(self) -> None:
        for relative in REQUIRED_DOCUMENTS:
            self.assertTrue((ROOT / relative).is_file(), relative)

    def test_documented_cli_examples_use_populated_csvs(self) -> None:
        # A copy-paste README example must actually run. Header-only templates
        # document the column contract but make the command fail with
        # "no compound rows found" — use shipped populated datasets instead.
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        inputs = re.findall(r"python -m pnn(?:\.rank)?\s+([\w./-]+\.csv)", readme)
        self.assertTrue(inputs, "no documented pnn CLI examples found")
        empty = []
        for rel in inputs:
            path = ROOT / rel
            if not path.is_file():
                continue  # illustrative placeholder (e.g. compounds.csv)
            lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if len(lines) < 2:
                empty.append(rel)
        self.assertEqual(empty, [], f"documented examples use header-only CSVs: {empty}")

    def test_python_version_requirement_is_consistent(self) -> None:
        # README is the gating onboarding doc; every other doc that states a
        # Python floor must agree with it or agents cold-start on the wrong version.
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        floor = re.search(r"Python (3\.\d+)\+", readme)
        self.assertIsNotNone(floor, "README.md has no 'Python 3.x+' prerequisite")
        for relative in REQUIRED_DOCUMENTS:
            if relative == "README.md":
                continue
            text = (ROOT / relative).read_text(encoding="utf-8")
            for other in re.findall(r"Python (3\.\d+)\+", text):
                self.assertEqual(
                    other, floor.group(1),
                    f"{relative} states Python {other}+ but README says {floor.group(1)}+",
                )

    def test_readme_tree_lists_only_existing_paths(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        heading = re.search(r"^#+ .*What's Inside", readme, re.M)
        self.assertIsNotNone(heading, "README.md has no What's Inside section")
        fence = readme.find("```text\n", heading.end())
        self.assertGreaterEqual(fence, 0, "README.md has no text tree")
        start = fence + len("```text\n")
        block = readme[start:readme.find("```", start)]
        stack = []
        checked = 0
        for line in block.splitlines():
            match = re.match(r"^((?:│   |    )*)(?:├── |└── )(.+?)(?:\s{2,}<-.*)?$", line)
            if match is None:
                continue
            depth = len(match.group(1)) // 4
            name = match.group(2).strip().rstrip("/")
            stack = stack[:depth] + [name]
            self.assertTrue((ROOT / "/".join(stack)).exists(), "/".join(stack))
            checked += 1
        self.assertGreater(checked, 0, "README.md tree lists no paths")

    def test_readme_tree_helper_matches_readme_tree(self) -> None:
        # update_readme_tree.py overwrites the README tree from NEW_TREE. If the
        # helper's template drifts, running it drops entries and breaks the docs
        # gate. The template must stay byte-for-byte in sync with README.
        from tools.update_readme_tree import NEW_TREE

        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        heading = re.search(r"^#+ .*What's Inside", readme, re.M)
        fence = readme.find("```text\n", heading.end())
        start = fence + len("```text\n")
        end = readme.find("```", start)
        actual = [line.rstrip() for line in readme[start:end].splitlines() if line.strip()]
        template = [line.rstrip() for line in NEW_TREE[len("```text\n") : NEW_TREE.rfind("```")].splitlines() if line.strip()]
        self.assertEqual(template, actual, "NEW_TREE drifted from the README tree")


if __name__ == "__main__":
    unittest.main()
