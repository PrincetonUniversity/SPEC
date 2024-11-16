import yaml

# Helper script to keep spec_refs (used to generate doxygen bibliography) 
# in sync with the machine readable CITATION.cff file. 
def cff_to_bibtex(cff_file, output_file):
    with open(cff_file, "r") as file:
        cff_data = yaml.safe_load(file)
    
    bibtex_entries = []

    for ref in cff_data.get("references", []):
        if ref.get("type") == "article":
            authors = " and ".join(
                f"{author['given-names']} {author['family-names']}" 
                for author in ref.get("authors", [])
            )
            bibtex_entry = (
                f"@article{{{ref.get('doi', 'unknown').replace('/', '_')},\n"
                f"  title={{ {ref.get('title')} }},\n"
                f"  author={{ {authors} }},\n"
                f"  journal={{ {ref.get('journal')} }},\n"
                f"  volume={{ {ref.get('volume')} }},\n"
                f"  number={{ {ref.get('issue', '')} }},\n"
                f"  pages={{ {ref.get('pages', '')} }},\n"
                f"  year={{ {ref.get('year')} }},\n"
                f"  doi={{ {ref.get('doi')} }},\n"
                f"  url={{ {ref.get('url')} }}\n"
                f"}}"
            )
            bibtex_entries.append(bibtex_entry)
    
    with open(output_file, "w") as file:
        file.write("\n\n".join(bibtex_entries))
    
    print(f"BibTeX entries written to {output_file}")

# Usage
cff_to_bibtex("CITATION.cff", "spec_refs.bib")
