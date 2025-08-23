import os

def main():
    # Pobierz ścieżkę pliku txt od użytkownika
    input_txt = input("Podaj ścieżkę do pliku txt z listą peptydów: ").strip()

    if not os.path.exists(input_txt):
        print(f"Błąd: plik {input_txt} nie istnieje.")
        return

    # Wczytaj peptydy z pliku txt
    with open(input_txt, "r", encoding="utf-8") as f:
        peptides = [line.strip() for line in f if line.strip()]

    if not peptides:
        print("Brak peptydów w pliku wejściowym.")
        return

    # Przygotuj folder wyjściowy
    output_dir = "/home/marta/Desktop/OPENBABEL_INPUT"
    os.makedirs(output_dir, exist_ok=True)

    # Wyciągnij nazwę pliku wejściowego (bez rozszerzenia)
    base_name = os.path.splitext(os.path.basename(input_txt))[0]

    # Utwórz nazwę pliku wynikowego (FASTA)
    output_fasta = os.path.join(output_dir, f"{base_name}.fasta")

    # Zapisz sekwencje w formacie FASTA
    with open(output_fasta, "w", encoding="utf-8") as f:
        for i, pep in enumerate(peptides, start=1):
            f.write(f">{base_name}_pep{i}\n")
            f.write(pep + "\n")

    print(f"Zapisano plik FASTA z {len(peptides)} peptydami: {output_fasta}")

if __name__ == "__main__":
    main()
