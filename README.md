# IEEE-WCL-LOCEDA
(c) 2024 Gangle Sun email: sungangleseu@gmail.com and sungangle@seu.edu.cn. 

This code implements the LOw-Coherence sEquence Design Algorithm (LOCEDA) from our WCL paper [1], which can generate sequences with (i) low-coherence, (ii) arbitrary lengths, (iii) any number of sequences, (iv) support for adaptable subcarrier assignments in orthogonal frequency-division multiple access (OFDMA) systems, and (v) compliance with user-defined PAPR constraints.

[1] G. Sun, W. Wang, W. Xu, and C. Studer, "Low-Coherence Sequence  Design Under PAPR Constraints," IEEE Wireless Commun. Lett., 2024.

This paper is also available at https://arxiv.org/abs/2407.21400v2

If you find our code and paper helpful, we would greatly appreciate it if you could cite our work. Thank you very much! ^o^

## IEEE BibTeX helper

Use `convert_bib.py` to convert BibTeX entries downloaded from IEEE Xplore into the abbreviated IEEE format (adds the month, wraps acronyms, and replaces journal names with the strings from `IEEEabrv.bib`).

```bash
python convert_bib.py input.bib --output output.bib
```

The script reads one or more entries from `input.bib` and writes the formatted result to stdout or to `output.bib` if provided.
