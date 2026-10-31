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


if __name__ == "__main__":
    unittest.main()
