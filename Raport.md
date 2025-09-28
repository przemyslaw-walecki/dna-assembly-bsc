# AI w Asemblacji de Novo - rozwiązywanie bąbelków
## Podsumowanie przeprowadzonych prac/eksperymentów

### 1. Podstawowy assembler
- asemblacja w Pythonie, przepisana na Rust (basic_dna_assembler -> rust_basic_dna_assembler)
- potwierdzone działanie na losowym genomie z deterministycznym sekwencjonowaniem 
- sprawdzone działanie na przykładowym genomie (drożdże), sekwenator Pirs
- podstawowe mechanizmy naprawy grafu: usuwanie martwych gałęzi & detekcja i rozwiązywanie bąbelków (metoda z pokryciem)
- ocena asemblacji przy pomocy minimap2 (alignment, avg_length)

### 2. Generowanie danych
- ekstrakcja grafu de Bruijn do pliku GFA, wykrywanie bąbelków za pomocą narzędzia BubbleGun
- etykietowanie bąbelków poprzez odnajdywanie ścieżki bąbelka w genomie referencyjnym, Reference-anchored walk, Minimal Unique Extension
- augmentacja danych - zmiana wartości pokrycia, zmiana losowego nukleotydu w k-merze

### 3. Użycie AI
- testowany model GNN + Feed Forward, dane podawane jako graf, wierzchołki = k-mery, połączenia = pokrycie
- zbierane wszystkie decyzje w bąbelku
- testowane na genomach bakterii, 1 genom tylko testowy (bąbelki nie były używane w treningu)
- najlepsze acc = 70% vs 50% heurystyka pokrycia
- oceniane tylko bąbelki poniżej określonego stostunku pokrycia (najlepsze wyniki dla poziomu <3.5)
- testowane różne hiperparametry modelu, rozmiary danych, siły augmentacji

### 4. Inne różne używane narzędzia
- Repozytorium kodu - Github
- Zarządzanie projektem - gałęzie, Issues
- Automatyczna dokumentacja - cargo docs + mkdocs
- proste skrypty shellowe/Python 
- Assembler w Python - poetry