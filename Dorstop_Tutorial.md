# Doorstop – Zarządzanie Wymaganiami Projektowymi

Niniejszy dokument przedstawia kompletny proces tworzenia i zarządzania wymaganiami projektowymi z wykorzystaniem narzędzia **Doorstop**. Doorstop to narzędzie typu open source, które umożliwia definiowanie, organizowanie i śledzenie wymagań w formacie tekstowym (YAML/Markdown), co czyni je doskonałym wyborem dla zespołów projektowych pracujących w środowiskach kontrolowanych przez systemy kontroli wersji (np. Git).

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

W wersji Doorstop v3.0 plik `.doorstop.yml` **nie jest tworzony automatycznie**.

Bez tego pliku Doorstop korzysta z wbudowanych ustawień, które zakładają:

- katalog z dokumentami: `docs/`
- rozszerzenie plików: `yml`
- brak dodatkowych ustawień niestandardowych

Aby uzyskać pełną kontrolę nad strukturą i sposobem działania narzędzia, należy utworzyć pli `.doorstop.yml`.

Przykładowa struktura:

```yaml
repo: docs
ext: yml
paths:
  - SYS
  - SW
encoding: utf-8
offline: true
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
| `ref`        | ciąg znaków lub lista | _(opcjonalne)_ Referencje zewnętrzne, np. normy, standardy (np. `ISO26262`, `RFC7519`), powiązane GitHub Issues |
| `derived`    | logiczny       | _(opcjonalne)_ Czy wymaganie zostało pochodnie wygenerowane z innego (`true`/`false`).  |
| `rationale`  | ciąg znaków    | _(opcjonalne)_ Uzasadnienie potrzeby danego wymagania.                                   |
| `parent`     | ciąg znaków    | _(przestarzałe)_ Stare pole używane przed `links`, obecnie niezalecane.                 |

### Dodatkowo, doorstop umożliwia podawanie dowolnych innych pól własnych, są one omijane przez silnik, ale mogą być używane przez inne narzędzia.

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
uid: SYS002
text: System powinien przechowywać dane w bezpieczny sposób.
level: 1
active: true
rationale: Ochrona danych użytkownika jest wymagana przez prawo oraz dobre praktyki bezpieczeństwa.
ref:
  - https://www.iso.org/isoiec-27001-information-security.html
  - GDPR-ART32
notes: Dotyczy zarówno danych osobowych, jak i danych systemowych przechowywanych lokalnie lub w chmurze.
file: docs/specs/data-security.md
test: untested
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

---
## Eksportowanie wymagań 
Doorstop udostępnia możliwość wyeksportowania wymagań do zbiorowych plików yaml, csv, tsv albo xlsx:
```bash
doorstop export <prefix|all> <plik> [-format]
```
prefix - nazwa dokumentu  
plik - ścieżka do pliku wyjściowego  
format - specyfikacja formatu (y/c/t/x)

---

### Pozostałe komendy

```bash
doorstop review <Nazwa> # Oznacza wymagania jako zrecenzowane przez <Nazwa>
doorstop unlink SW001 SYS001
```
---
## Przykłady użycia - integracja z CI/CD w GitHubie:
Doorstop bardzo dobrze sprawdza się jako dosyć wszechstronne narzędzie umożliwiające integrację z różnymi zewnętrznymi systemami. Poniżej przykład integracji z CI/CD w GitHubie za pomocą skryptu w Pythonie w wersji Doorstop v1.2:
Plik `check_requirements.py`
```python
import sys
from doorstop.core import tree

REQUIRED_PREFIX = "SYS"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}

def main():
    t = tree.build()

    doc = t.get_document(REQUIRED_PREFIX)
    if not doc:
        print(f"ERROR: Document '{REQUIRED_PREFIX}' not found.", file=sys.stderr)
        sys.exit(1)

    active_items = [item for item in doc.items if item.active]
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
          python-version: 3.10
      - name: Install dependencies
        run: pip install doorstop
      - name: Check requirement coverage
        run: python check_requirements_tested.py

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
Możliwe było wystawienie lokalnego API do zarządzania wymaganiami przez zapytania REST. Funkcjonalność została usunięta w v3.0.
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