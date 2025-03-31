import sys
import yaml

EXPORT_FILE = "requirements.yaml"
REQUIRED_PREFIX = "SYS"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}

def main():
    with open(EXPORT_FILE, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if not data:
        print("ERROR: Exported YAML is empty or invalid.", file=sys.stderr)
        sys.exit(1)

    relevant_items = [
        item for uid, item in data.items()
        if uid.startswith(REQUIRED_PREFIX) and item.get("active", True)
    ]

    if not relevant_items:
        print(f"ERROR: No active requirements with prefix '{REQUIRED_PREFIX}'.", file=sys.stderr)
        sys.exit(1)

    tested = sum(
        1 for item in relevant_items
        if isinstance(item.get("test"), str) and item["test"].lower() in VALID_TEST_STATUSES
    )

    total = len(relevant_items)
    coverage = tested / total

    print(f"{REQUIRED_PREFIX} coverage: {tested}/{total} ({coverage:.1%})")

    if coverage < REQUIRED_COVERAGE:
        print(f"ERROR: Test coverage below threshold ({REQUIRED_COVERAGE:.0%}).", file=sys.stderr)
        sys.exit(1)

    print("SUCCESS: Requirement test coverage OK.")
    sys.exit(0)

if __name__ == "__main__":
    main()
