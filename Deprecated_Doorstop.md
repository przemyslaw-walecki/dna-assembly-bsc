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