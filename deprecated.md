## Tabela funkcjonalności

| Funkcja                  | Komenda lub plik                      |
|--------------------------|----------------------------------------|
| Inicjalizacja            | `doorstop create`                     |
| Dodanie dokumentu        | `doorstop create SYS`                 |
| Dodanie wymagania        | `doorstop add SYS`                    |
| Tworzenie powiązań       | `doorstop create SW -p SYS`           |
| Drzewo wymagań           | `doorstop tree`                       |
| Walidacja                | `doorstop verify`                     |
| Eksport                  | `doorstop export --format markdown`   |
| Szukanie                 | `doorstop find {--pole} {wartosc1} {wartosc2}`                 |
| Status testów            | `test: implemented` w pliku YAML      |

## Eksport dokumentacji

Doorstop umożliwia eksport wymagań do różnych formatów:

```bash
doorstop export --format markdown > wymagania.md
doorstop export --format json > wymagania.json
doorstop export --format html > wymagania.html
```

### Lista wymagań: `doorstop list`


Polecenie `doorstop list` pozwala na wyświetlenie wszystkich wymagań zdefiniowanych w projekcie. Dodatkowo możliwe jest filtrowanie według wybranych atrybutów:
```bash
doorstop list {--pole} {wartość}
doorstop list --test untested
doorstop list --{active/inactive} # Tylko wymagania z danym statusem
doorstop list --doc SW # Tylko dany dokument
```
## Przeszukiwanie i filtrowanie wymagań

Narzędzie Doorstop umożliwia przeszukiwanie oraz filtrowanie dokumentów z wymaganiami za pomocą wbudowanych poleceń CLI.


Polecenie `doorstop find` pozwala wyszukiwać wyrażenia w polu `text` lub `notes`
```bash 
doorstop find login authentication # Wiele słów kluczowych po spacji
```
