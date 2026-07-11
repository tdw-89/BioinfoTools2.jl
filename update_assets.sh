#!/usr/bin/env sh
curl -o ./assets/SOFA.json https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/heads/master/Ontology_Files/subsets/SOFA.json
echo "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/heads/master/Ontology_Files/subsets/SOFA.json" > ./assets/SOFA.json.info
echo $(date) >> ./assets/SOFA.json.info