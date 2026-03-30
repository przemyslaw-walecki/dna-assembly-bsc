# dna-assembly-bsc

## Instrukcja użytkownika

### 1. Uruchomienie przepływu - od odczytu do decyzji
Aby uruchomić pełny przepływ asemblujący genom zadanej bakterii w oparciu o predefiniowany model GNN, należy uruchomić skrypt run_pipeline.sh  
  
Wymagania:
- Odczyty genomu w formacie fastq wygenerowane przy użyciu sekwenatora pirs
- Python w wersji > 3.10

### 2. Uruchomienie eksperymentów - folder bubble_resolution_gnn
Aby uruchomić poszczególne eskperymenty, należy uruchomić poszczególne pliki .ipynb znajdujące się w folderze bubble_resolution_gnn.
  
Wymagania:
- Python > 3.10 
- torch, torch_geometric

Zalecane:  
- cuda
