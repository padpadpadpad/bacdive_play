PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX d3o: <https://purl.dsmz.de/schema/>
  PREFIX bd: <https://purl.dsmz.de/bacdive/>
  SELECT DISTINCT ?bacdiveid ?length
WHERE {
  ?morphology_or_size_object d3o:hasLength ?length . # Find the object that has a length
  
  # This morphology/size object describes a specific strain
  ?morphology_or_size_object d3o:describesStrain ?strain .
  
  # The strain then has the BacDive ID
  ?strain d3o:hasBacDiveID ?bacdiveid .
}