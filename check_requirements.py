import sys
from doorstop.core.project import load

REQUIRED_PREFIX = "SYS"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}

def main():
    project = load()

    sys_doc = next((d for d in project.documents if d.prefix == REQUIRED_PREFIX), None)
    if not sys_doc:
        print(f"ERROR: Document '{REQUIRED_PREFIX}' not found.", file=sys.stderr)
        sys.exit(1)

    active_items = [item for item in sys_doc.items if item.active]
    if not active_items:
        print(f"ERROR: No active items in document '{REQUIRED_PREFIX}'.", file=sys.stderr)
        sys.exit(1)

    tested_items = [
        item for item in active_items
        if isinstance(item.test, str) and item.test.lower() in VALID_TEST_STATUSES
    ]

    total = len(active_items)
    tested = len(tested_items)
    coverage = tested / total

    print(f"{REQUIRED_PREFIX} coverage: {tested}/{total} ({coverage:.1%})")

    if coverage < REQUIRED_COVERAGE:
        print(f"ERROR: Test coverage below threshold ({REQUIRED_COVERAGE:.0%}).", file=sys.stderr)
        sys.exit(1)

    print("SUCCESS: Requirement test coverage OK.")
    sys.exit(0)

if __name__ == "__main__":
    main()
