import pandas as pd
import os

def main():
    # Pobierz ścieżkę pliku txt od użytkownika
    input_txt = input("Podaj ścieżkę do pliku txt z listą peptydów: ").strip()

    if not os.path.exists(input_txt):
        print(f"Błąd: plik {input_txt} nie istnieje.")
        return

    # Ścieżka do pliku .ods z referencyjnymi peptydami
    reference_ods = "/home/marta/Desktop/known_ace_inhibitory_peptides.ods"

    if not os.path.exists(reference_ods):
        print(f"Błąd: plik {reference_ods} nie istnieje.")
        return

    # Wczytaj znane peptydy z pliku .ods (pierwsza kolumna, pomijamy header)
    known_df = pd.read_excel(reference_ods, engine="odf")
    known_peptides = set(known_df.iloc[:, 0].dropna().astype(str).str.strip())

    # Wczytaj peptydy z pliku txt
    with open(input_txt, "r", encoding="utf-8") as f:
        input_peptides = {line.strip() for line in f if line.strip()}

    # Znajdź nowe peptydy
    new_peptides = input_peptides - known_peptides

    # Przygotuj folder wyjściowy
    output_dir = "/home/marta/Desktop/novel_peptides_compare"
    os.makedirs(output_dir, exist_ok=True)

    # Wyciągnij nazwę gatunku z pliku wejściowego (bez rozszerzenia)
    species = os.path.splitext(os.path.basename(input_txt))[0]

    # Utwórz nazwę pliku wynikowego
    output_file = os.path.join(output_dir, f"{species}_potentially_novel_peptides.txt")

    # Zapisz wynik do nowego pliku
    with open(output_file, "w", encoding="utf-8") as f:
        for pep in sorted(new_peptides):
            f.write(pep + "\n")

    print(f"Zapisano {len(new_peptides)} nowych peptydów do pliku: {output_file}")

if __name__ == "__main__":
    main()
