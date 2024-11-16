import yaml

# Helper script to keep spec_refs (used to generate doxygen bibliography) 
# in sync with the machine readable CITATION.cff file. 
def cff_to_bibtex(cff_file, output_file):
    with open(cff_file, "r") as file:
        cff_data = yaml.safe_load(file)
    
    bibtex_entries = []
    used_keys = set()

    for ref in cff_data.get("references", []):
        if ref.get("type") == "article":
            # Extract year and first author's last name
            year = ref.get("year", "unknown")
            first_author_last_name = (
                ref.get("authors", [{}])[0].get("family-names", "unknown").lower()
            )
            
            # Construct the initial BibTeX key
            base_key = f"y{year}_{first_author_last_name}"
            key = base_key
            
            # Ensure the key is unique
            count = 1
            while key in used_keys:
                key = f"{base_key}_{count}"
                count += 1
            used_keys.add(key)
            
            # Construct author string
            authors = " and ".join(
                f"{author['given-names']} {author['family-names']}" 
                for author in ref.get("authors", [])
            )
            
            bibtex_entry = (
                f"@article{{{key},\n"
                f"  title={{ {ref.get('title')} }},\n"
                f"  author={{ {authors} }},\n"
                f"  journal={{ {ref.get('journal')} }},\n"
                f"  volume={{ {ref.get('volume')} }},\n"
                f"  number={{ {ref.get('issue', '')} }},\n"
                f"  pages={{ {ref.get('pages', '')} }},\n"
                f"  year={{ {year} }},\n"
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
