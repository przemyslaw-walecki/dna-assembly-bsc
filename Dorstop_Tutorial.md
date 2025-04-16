# Doorstop – Zarządzanie Wymaganiami Projektowymi

Doorstop to narzędzie typu open source, które umożliwia definiowanie, organizowanie i śledzenie wymagań w formacie tekstowym YAML.
---

## Instalacja Doorstop

```bash
pip install doorstop
```

---

## Inicjalizacja projektu

Podstawowa struktura projektu:

```bash
mkdir doorstop-tutorial
cd doorstop-tutorial
git init
mkdir docs
```

---
## Plik konfiguracyjny `.doorstop.yml`

W wersji v3.0 Doorstop korzysta z wbudowanych ustawień, które zakładają:

- katalog z dokumentami: `docs/`
- rozszerzenie plików: `yml`
- brak dodatkowych ustawień niestandardowych

Aby uzyskać pełną kontrolę nad strukturą i sposobem działania narzędzia, należy edytować plik `.doorstop.yml` wewnątrz folderu dokumentu.

Przykładowa struktura:

```yaml
settings:
  digits: 3
  itemformat: yaml
  parent: WT
  prefix: WF
  sep: ''
```
---

## Struktura wymagań: podział na WF i WT

W ramach poniższego dokumentu zastosowano dwupoziomowy podział:

| Typ dokumentu | Opis                                                                 |
|---------------|----------------------------------------------------------------------|
| **WF**        | Wymagania funkcjonalne. Opisują, co system powinien robić.           |
| **WT**        | Wymagania testowe. Określają sposób weryfikacji spełnienia wymagań. |

Doorstop umożliwia reprezentowanie wymagań jako elementów powiązanych ze sobą poprzez relacje nadrzędny–podrzędny lub wiele-do-wielu (graf wymagań).


## Tworzenie dokumentu z wymaganiami funkcjonalnymi (WF)

Wymagania funkcjonalne reprezentują poszczególne funkcjonalości systemu.

```bash
doorstop create WF docs/WF # doorstop create {nazwa} {folder_na_wymagania}
doorstop add WF # Dostępne dodatkowe flagi: [-c ilość wymagań] [--file plik_wejściowy] [--template własny_template]  
```
Każde wymaganie jest reprezentowane jako plik YAML. Poniżej opisano możliwe pola w pliku wymagań:

| Pole         | Typ           | Opis                                                                                     |
|--------------|----------------|------------------------------------------------------------------------------------------|
| `uid`        | ciąg znaków    | Unikalny identyfikator wymogu (np. `WF001`, `WT003`).                                  |
| `text`       | ciąg znaków    | Treść wymagania w języku naturalnym – opisuje oczekiwane zachowanie/system.             |
| `level`      | liczba         | Poziom wymagania wewnątrz grupy: 1.0 2.0 - ten sam poziom, 1.1 - pod 1.0, widoczne poprzez wcięcia po zpublikowaniu do html         |
| `active`     | logiczny       | Czy wymaganie jest aktywne (`true`) czy nieaktywne (`false`).                            |
| `links`      | lista          | _(opcjonalne)_ Lista identyfikatorów powiązanych wymagań nadrzędnych (np. `WF001`).    |
| `test`       | ciąg znaków    | _(opcjonalne)_ Status testu – np. `implemented`, `untested`, `passed`, `failed`.        |
| `file`       | ciąg znaków    | _(opcjonalne)_ Ścieżka do pliku powiązanego z wymaganiem (np. kod, test, dokument).     |
| `notes`      | ciąg znaków    | _(opcjonalne)_ Dodatkowe uwagi, wyjaśnienia, wymagania nieformalnie opisane.            |
| `ref`        | ciąg znaków lub lista | _(opcjonalne)_ Referencje zewnętrzne, z innych plików w projekcie, Doorstop sam znajdzie ich lokacje |
| `derived`    | logiczny       | _(opcjonalne)_ Czy wymaganie zostało pochodnie wygenerowane z innego (`true`/`false`).  |
| `rationale`  | ciąg znaków    | _(opcjonalne)_ Uzasadnienie potrzeby danego wymagania.                                   |
| `parent`     | ciąg znaków    | _(przestarzałe)_ Stare pole używane przed `links`, obecnie niezalecane.                 |

Dodatkowo, doorstop umożliwia podawanie dowolnych innych pól własnych, są one omijane przez silnik, ale mogą być używane przez własne narzędzia.

### Doorstop umożliwia podstawowe zarządzanie wymaganiami:
```bash 
doorstop edit WF001 # Otwiera plik yaml w domyślnym edytorze
doorstop remove WF002
doorstop reorder WF # Uporządkowanie UID, np po remove
```
---

Przykładowy plik `docs/WF/WF001.yml`:

```yaml
active: true
derived: false
header: ''
level: 1.0
links: []
normative: true
ref: Aut
reviewed: uMYtrrGresKU4sLxHaHzH6n-vj7K1_pcAX7fPpIfqGQ=
text: |
  System powinien umożliwić logowanie użytkownika.
uid: WF001
```

Dodanie kolejnego wymagania:

```bash
doorstop add WF
```
W polu ref można umieszczać tekst. Doorstop automatycznie przeszuka katalog roboczy i odszuka wpisany tekst w plikach. Widoczne po np zpublikowaniu wymagań do formatu html.

---

## Łączenie grup wymagań

Doorstop umożliwia tworzenie powiązań między wymaganiami, wskazując ich zależność lub pochodzenie.


```bash
doorstop create WT docs/WT
doorstop add WT
doorstop link WT001 WF001 # Powiązanie widoczne tylko w WT001
```

Przykładowy plik `WT001.yml`:

```yaml
uid: WT001
text: Użytkownik otrzymuje klucz sesji po podaniu poprawnych danych logowania
level: 1.0
links:
  - WF001: null
active: true
test: implemented
```
## Walidacja struktury wymagań

Aby sprawdzić poprawność struktury naszych wymagań, należy po prostu uruchomić:
```bash
doorstop
```
Wówczas, jeżeli struktura wymagań jest poprawna, Doorstop "podpisze" odpowiednie pola w wymaganiu hashem SHA-256

---
## Eksportowanie wymagań 
Doorstop udostępnia możliwość wyeksportowania wymagań do zbiorowych plików yaml, csv, tsv albo xlsx:
```bash
doorstop export <prefix|all> <plik> [-format]
```
prefix - nazwa dokumentu  
plik - ścieżka do pliku wyjściowego  
format - specyfikacja formatu (y/c/t/x)

Można również publikować dokumenty jako pliki html:
```bash
doorstop publish <prefix|all> [ścieżka]
```
---

### Pozostałe komendy

```bash
doorstop review <prefix|all> # Oznacza wymagania jako zrecenzowane - zostaje ustawiony hash SHA-256 na podstawie treści wymagania
doorstop unlink WT001 WF001
```
---
## Przykłady użycia - integracja z CI/CD w GitHubie:
Doorstop bardzo dobrze sprawdza się jako dosyć wszechstronne narzędzie umożliwiające integrację z różnymi zewnętrznymi Systemami. Poniżej przykład integracji z CI/CD w GitHubie za pomocą skryptu w Pythonie:
Plik `check_requirements.py`
```python
import WF
import yaml

EXPORT_FILE = "requirements.yaml"
REQUIRED_PREFIX = "WF"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}

def main():
    with open(EXPORT_FILE, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if not data:
        print("ERROR: Exported YAML is empty or invalid.", file=WF.stderr)
        WF.exit(1)

    relevant_items = [
        item for uid, item in data.items()
        if uid.startWTith(REQUIRED_PREFIX) and item.get("active", True)
    ]

    if not relevant_items:
        print(f"ERROR: No active requirements with prefix '{REQUIRED_PREFIX}'.", file=WF.stderr)
        WF.exit(1)

    tested = sum(
        1 for item in relevant_items
        if isinstance(item.get("test"), str) and item["test"].lower() in VALID_TEST_STATUSES
    )

    total = len(relevant_items)
    coverage = tested / total

    print(f"{REQUIRED_PREFIX} coverage: {tested}/{total} ({coverage:.1%})")

    if coverage < REQUIRED_COVERAGE:
        print(f"ERROR: Test coverage below threshold ({REQUIRED_COVERAGE:.0%}).", file=WF.stderr)
        WF.exit(1)

    print("SUCCESS: Requirement test coverage OK.")
    WF.exit(0)

if __name__ == "__main__":
    main()

```
Plik `.github/workflows/verify_requirements.yml`
```yaml
name: Check Requirement Coverage

on: [push, pull_request]

jobs:
  test-coverage:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.10'

      - name: Install dependencies
        run: pip install doorstop pyyaml

      - name: Export requirements
        run: doorstop export -y WF requirements.yaml


      - name: Check requirement coverage
        run: python check_requirements.py

```

---

### Źródła
https://github.com/doorstop-dev/doorstop
https://doorstop.readthedocs.io/en/latest/