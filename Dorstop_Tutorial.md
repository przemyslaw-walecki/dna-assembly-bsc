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
repo: docs
ext: yml
paths:
  - SYS
  - SW
encoding: utf-8
verbose: true
index_name: index
```
---

## Tworzenie dokumentu z wymaganiami systemowymi (SYS)

Wymagania systemowe reprezentują ogólne cele i założenia projektu.

```bash
doorstop create SYS docs/SYS # doorstop create {nazwa} {folder_na_wymagania}
doorstop add SYS # Dostępne dodatkowe flagi: [-c ilość wymagań] [--file plik_wejściowy] [--template własny_template]  
```

Każde wymaganie jest reprezentowane jako plik YAML. Poniżej opisano możliwe pola w pliku wymagań:

| Pole         | Typ           | Opis                                                                                     |
|--------------|----------------|------------------------------------------------------------------------------------------|
| `uid`        | ciąg znaków    | Unikalny identyfikator wymogu (np. `SYS001`, `SW003`).                                  |
| `text`       | ciąg znaków    | Treść wymagania w języku naturalnym – opisuje oczekiwane zachowanie/system.             |
| `level`      | liczba         | Poziom wymagania – 1 = ogólne, wysokopoziomowe, 2+ = coraz bardziej szczegółowe.         |
| `active`     | logiczny       | Czy wymaganie jest aktywne (`true`) czy nieaktywne (`false`).                            |
| `links`      | lista          | _(opcjonalne)_ Lista identyfikatorów powiązanych wymagań nadrzędnych (np. `SYS001`).    |
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
doorstop edit SYS001 # Otwiera plik yaml w domyślnym edytorze
doorstop remove SYS002
doorstop reorder SW # Uporządkowanie UID, np po remove
```
---

Przykładowy plik `docs/SYS/SYS001.yml`:

```yaml
uid: SYS001
text: System powinien obsługiwać uwierzytelnianie użytkownika.
level: 1
active: true
```

Dodanie kolejnego wymagania:

```bash
doorstop add SYS
```

Przykładowa konfiguracja `docs/SYS/SYS002.yml`:

```yaml
active: true
derived: false
file: docs/specs/data-security.md
header: |
  To jest testowy Header
level: 1.1
links: []
normative: true
notes: |
  'Dotyczy zarówno danych osobowych, jak i danych systemowych przechowywanych lokalnie lub w chmurze.'
rationale: |
  'Ochrona danych użytkownika jest wymagana przez prawo oraz dobre praktyki bezpieczeństwa.'
ref: www.iso.org/isoiec-27001-information-security.html
reviewed: n-1OZqmSlOe8UqIrKcueY-zYw-SFDmifkr9o1IS676A=
test: ok
text: |
  System powinien przechowywać dane w bezpieczny sposób.
uid: SYS002
```

---

## Łączenie grup wymagań

Doorstop umożliwia tworzenie powiązań między wymaganiami, wskazując ich zależność lub pochodzenie.


```bash
doorstop create SW docs/SW
doorstop add SW
doorstop link SW001 SYS001 # Powiązanie widoczne tylko w SW001
```

Przykładowy plik `SW001.yml`:

```yaml
uid: SW001
text: Aplikacja powinna implementować formularz logowania.
level: 2
links:
  - SYS001: null # Możliwe ustawienie dodatkowych atrbytuów połączenia, np. typ lub uzasadnienie
active: true
test: implemented
```
## Walidacja struktury wymagań

Aby sprawdzić poprawność struktury naszych wymagań, należy po prostu uruchomić:
```bash
doorstop
```

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
doorstop review <prefix|all> # Oznacza wymagania jako zrecenzowane - zostaje ustawiony unikatowy klucz SHA-256 na podstawie treści wymagania
doorstop unlink SW001 SYS001
```
---
## Przykłady użycia - integracja z CI/CD w GitHubie:
Doorstop bardzo dobrze sprawdza się jako dosyć wszechstronne narzędzie umożliwiające integrację z różnymi zewnętrznymi systemami. Poniżej przykład integracji z CI/CD w GitHubie za pomocą skryptu w Pythonie:
Plik `check_requirements.py`
```python
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
        run: doorstop export -y SYS requirements.yaml


      - name: Check requirement coverage
        run: python check_requirements.py

```

---

## Poprzednie wersje:

Doorstop v1.x udostępniał znacznie więcej funkcji do zarządzania i nadzorowania wymagań.  
Większość z nich zaczęła być uznawana za przestarzała i przesuwana do zewnętrznych pluginów w wersji v2.0.

Przykładowe komendy dostępne w v1.x, w tym różnice implementacji:

| Funkcja                  | Komenda lub plik                      |
|--------------------------|----------------------------------------|
| Inicjalizacja            | `doorstop create`                     |
| Dodanie dokumentu        | `doorstop create SYS`                 |
| Dodanie wymagania        | `doorstop add SYS`                    |
| Tworzenie powiązań       | `doorstop create SW -p SYS`           |
| Drzewo wymagań           | `doorstop tree`                       |
| Walidacja struktury wymagań           | `doorstop verify`       |
| Eksport                  | `doorstop export --format markdown`   |
| Szukanie                 | `doorstop find {--pole} {wartosc1} {wartosc2}`                 |
| Filtrowanie              | `doorstop list {--pole} {wartość}` 
|                           | `doorstop list --{active/inactive}` - status  |

## Serwer Doorstop
Możliwe było wystawienie lokalnego API do zarządzania wymaganiami przez zapytania REST. Funkcjonalność raczej eksperymentalna, brak wsparcia.
```bash 
doorstop-server # Domyślnie uruchamia na localhost:7867
```


| Metoda   | Endpoint                         | Opis                                                                 |
|----------|----------------------------------|----------------------------------------------------------------------|
| `GET`    | `/docs`                          | Zwraca listę wszystkich dokumentów (np. `SYS`, `SW`).               |
| `GET`    | `/docs/<prefix>`                 | Zwraca metadane dokumentu (np. `index.yml` dla `SYS`).              |
| `GET`    | `/items`                         | Zwraca listę wszystkich wymagań ze wszystkich dokumentów.           |
| `GET`    | `/items/<uid>`                   | Zwraca szczegóły konkretnego wymagania.                             |
| `POST`   | `/items/<uid>`                   | Tworzy lub aktualizuje wymaganie (JSON w `body`).                   |
| `DELETE` | `/items/<uid>`                   | Usuwa wymaganie (jeśli obsługiwane przez wersję serwera).           |
| `GET`    | `/tree`                          | Zwraca strukturę drzewa powiązań między dokumentami.                |
| `GET`    | `/verify`                        | Weryfikuje poprawność powiązań i kompletność wymagań (jeśli dostępne). |
| `GET`    | `/export/markdown`               | Eksportuje dokumentację w formacie Markdown.                        |
| `GET`    | `/export/html`                   | Eksportuje dokumentację w formacie HTML (jeśli obsługiwane).        |
| `GET`    | `/export/pdf`                    | Eksportuje dokumentację do PDF (jeśli zintegrowano Pandoc/LaTeX).   |

### Źródła
https://github.com/doorstop-dev/doorstop
https://doorstop.readthedocs.io/en/latest/