# `gbk2tbl` – GenBank to NCBI TBL Converter

This tool converts **GenBank (.gbk)** files into **NCBI-style 5-column tab-delimited feature tables (.tbl)**, suitable for submission to GenBank via:

- [BankIt](https://www.ncbi.nlm.nih.gov/WebSub/?tool=bankit)
- [`tbl2asn`](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/)
- [GenBank Submission Portal](https://submit.ncbi.nlm.nih.gov/)

It supports:
- Proper **1-based coordinate conversion**
- Detection of **partial features** (`<`, `>`)
- Automatic inference of `codon_start` for **partial CDSs**
- Enhanced validation:
  - Gene-CDS consistency
  - Protein-coding gene detection
  - Location mismatch checks
- Optional qualifiers:
  - `transl_table`
  - `REFERENCE` (via PubMed ID)
  - `[offset=N]` directive for multi-part submissions

---

## 📦 Features

| Feature | Description |
|--------|-------------|
✅ **NCBI Compliance** | Fully conforms to NCBI feature table format guidelines |
✅ **Partial Features** | Supports `<start` and `>end` syntax |
✅ **Codon Start Inference** | Automatically computes `codon_start` for partial CDSs |
✅ **Gene-CDS Consistency Checks** | Detects mismatches between gene and CDS locations |
✅ **Validation & Logging** | Warnings for missing products, invalid genes, etc. |
✅ **Custom Qualifier Handling** | Sanitizes values, wraps long lines |
✅ **Reverse Strand Support** | Correctly handles features on the reverse strand |

---

## 🚀 Installation

### Clone the repo:
```bash
git clone https://github.com/MrGhost-Aymen/gbk2tbl.git
cd gbk2tbl
```

### Install dependencies:
```bash
pip install biopython
```

---

## 🛠️ Usage

```bash
python gbk_to_tbl.py input.gbk output.tbl [OPTIONS]
```

### 🔧 Options

| Option | Description |
|--------|-------------|
`--transl_table N` | Add transl_table qualifier (e.g., `2` for mitochondrial genomes) |
`--pmid PMID` | Add a REFERENCE line using a PubMed ID |
`--no-consistency-check` | Skip gene/CDS consistency checks |
`--offset N` | Add an offset to all positions after this file |
`-v, --verbose` | Enable debug-level logging |
`--version` | Show version info |
`--help` | Show help message |

### ✅ Example Commands

#### Basic usage:
```bash
python gbk_to_tbl.py input.gbk output.tbl
```

#### For mitochondrial genomes:
```bash
python gbk_to_tbl.py input.gbk output.tbl --transl_table 2
```

#### With literature reference:
```bash
python gbk_to_tbl.py input.gbk output.tbl --transl_table 2 --pmid 12345678
```

#### For multipart genome assembly:
```bash
python gbk_to_tbl.py part1.gbk part1.tbl --offset 0
python gbk_to_tbl.py part2.gbk part2.tbl --offset 50000
```

---

## 🧾 Output Format

The generated `.tbl` file follows NCBI’s **5-column tab-delimited feature table** format:

```
>Feature NC_012345
1	1000	gene
		gene	ATP6
1	1000	CDS
		product	ATP synthase subunit 6
		codon_start	1
```

Supports:
- Multi-exon features
- Trans-spliced genes
- Reverse-strand annotations
- Long qualifier wrapping

---

## 📁 Input Requirements

- Valid **GenBank (.gbk)** file(s)
- Must include sequence and annotation data
- Supported feature types:  
  `gene`, `CDS`, `mRNA`, `tRNA`, `rRNA`, `ncRNA`, `exon`, `intron`, `repeat_region`, and more

---

## ⚠️ Validation & Warnings

The script includes smart validation:
- Warns about **missing CDSs** for protein-coding genes (like `nad1`, `cox2`, etc.)
- Detects **gene-CDS location mismatches**
- Flags **missing product qualifiers** in CDS features
- Validates sequence quality (e.g., short sequences, high N-count)

Example warning:
```
WARNING: Gene-CDS location mismatch for 'nad1': gene(195200-195585) vs CDS(13824-195585). Gene should span entire CDS range.
```

---

## 🧪 Example

Input:
```bash
python gbk_to_tbl.py example.gbk example.tbl --transl_table 2 --pmid 34567890
```

Output:
```text
>Feature NC_012345
1	1000	gene
		gene	ATP6
1	1000	CDS
		gene	ATP6
		product	ATP synthase subunit 6
		codon_start	1
		transl_table	2
1	1000	tRNA
		gene	trnM-CAU
		product	tRNA-Met
1	1000	REFERENCE
		PubMed	34567890
```

---

## 🧩 Why Use This Tool?

- ✅ Ready for NCBI submission
- 🧹 Clean, well-documented code
- 🧪 Built-in validation and logging
- 📁 Designed for batch processing and pipelines
- 🌐 MIT License – free to use and modify

---

## 📄 NCBI Guidelines Reference

This tool was built following the official **NCBI Feature Table Guidelines**, ensuring full compatibility with BankIt, tbl2asn, and GenBank Submission Portal.

📄 Full NCBI Feature Table Guide included in project documentation.

---

## 🧑‍💻 Contributing

Contributions are welcome! If you'd like to improve or extend this tool, please open an issue or submit a PR.

---

## 📬 Questions or Feedback?

Feel free to open an issue or reach out if you need:
- Custom feature support
- Integration into your annotation pipeline
- Batch processing scripts
- GUI version

---

## 📣 Acknowledgments

Special thanks to the NCBI team for their comprehensive documentation and standards.

---

**Ready to convert your `.gbk` files into NCBI-ready `.tbl` files?**

👉 Just run:  
```bash
python gbk_to_tbl.py your_file.gbk your_output.tbl --transl_table 2 
```

🚀 Happy submitting to GenBank!
