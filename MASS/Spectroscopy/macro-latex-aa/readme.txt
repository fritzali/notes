%-----------------------------------------------------------------------
% AA LaTeX class for Astronomy & Astrophysics						 9.4
%
%                                                 (c) EDP Sciences, 2026
%                                            tex-support@edpsciences.org
%-----------------------------------------------------------------------

Documentation:
https://www.aanda.org/for-authors/latex-issues/texnical-background-information

The following are part of the AA macro package, compressed into a single
archive: macro-latex-aa.zip

- macros directory
- aa_example.zip: Example of an article (LaTeX source)
- aa_example.pdf: Example of an article (PDF output)
- readme.txt: this file


%---------------------------------------------------------------------
%														     Changelog 
%---------------------------------------------------------------------

9.4 March 2026
- Introduction of the [corrauth] command to automatically put
the corresponding author address as footnote
- [email] command is reserved to extraction and will be ignored
during the PDF compilation

9.3 October 2025
- Behaviour of the [longauth] command is modified to shift the whole
author list and affiliations after references.

9.2 September 2024
- Fix affilation overflow with clearpage (when using longauth)
- Fix hyperlinks
- Flat pagination 
- Remove line numbering style
- Remove legacy orcid macro
- Add grffile style to manage complex strings for figure names