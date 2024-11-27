# README pipeline BNGP

## Wat doet de pipeline?
Deze pipeline is gemaakt voor onderzoek naar de bacterie *Listeria monocytogenes*. De pipeline is getest op 43 Oxford Nanopore reads. Deze data zijn verkregen op het
LCAB voor een onderzoek naar de aanwezigheid van deze pathogene bacterie in de voedselverwerkende industrie. Om de herkomst van uitbraken in kaart te brengen, is het 
belangrijk om de verwantschap tussen de *L. Monocytogenes* samples te bepalen. Door deze resulaten valt te herleiden of de uitbraak door een stam wordt veroorzaakt of
door verschillende onafhankelijke stammen. In deze pipeline wordt de verwantschap op twee verschillende wijzen bepaald: op basis van SNP's, en op basis van whole genome
hash sketches. De output van deze pipeline bestaat uit twee fylogenetische bomen (een van SNP en de ander van de hash sketches). 

## PIPELINE UITVOEREN:
Om de pipeline van BNGP uit te voeren, gebruik je de volgende command:

`snakemake --cores <aantal beschikbare cores> --resources lock=1`

## CONFIGURATIE:
In de config.yaml kan je 6 waardes aanpassen:
- Bij DIRECTORY_FQ moet je het pad opgeven waar de raw reads zitten.
- Bij REFGENOME en REFGENOME_BASE_NAME moet je het pad naar het referentiegenoom opgeven en de basename van het referentiegenoom (dus zonder pad).
- Voor het bewerken en filteren van de reads, zijn er 3 parameters.
- TRIM_BASES bepaalt hoeveel basen van je reads aan beide kanten worden verwijderd.
- MIN_PHRED bepaalt welke gemiddelde PHRED waarde je reads minimaal moeten hebben.
- MIN_LENGTH bepaalt wat de minimale lengte van je reads moet zijn.

## OUTPUT FILES:
Uit de pipeline komen de volgende output files:
- qc_results/ -> De resultaten van de QC. In het mapje 'old' staan de QC's van de oude reads, in edited staan de QC's van de bewerkte reads.
- edited_fastq/ -> Mapje met de bewerkte reads.
- BAM/ -> De gesorteerde bam bestanden.
- (refgenome).fai -> fasta-index van het refgenome.
- pileup/ -> Alle pileup files.
- bcf/ -> Hierin zitten de variant-calls, de variant-calls gefilterd op SNPS, bcf variant van het gefiltere vcf bestand en het geindexeerde bcf bestand.
- consensus/ -> De consensus files per sample, de consesus files met een aangepaste header en combined.fasta bevat alle consensus files bij elkaar. Hiervoor is het belangrijk dat alle consensus files dezelfde lengte hebben.
- SNP.* -> Alle files die afkomen van iqtree2
- SNP.treefile -> De fylogenetische SNP boom.
- flye/(sample)/assembly.fasta -> De assemblies uitgevoerd met flye.
- hashing/ -> Hash sketches van alle flye assemblies.
- comparison.csv -> Distance matrix tussen alle hash sketches.
- flye_tree.newick -> De fylogenetische flye/de novo boom.

## UITLEG SCRIPTS:
### quality_control.py:
Dit script stelt een QC rapport op. Hij krijgt een fastq bestand als input en hij geeft een JSON bestand als output.
Als eerste wordt het fastq bestand ingelezen. Het script blijft per read lezen, totdat hij door het fastq bestand heen is. De volgende checks worden in chronologische volgorde berekent.
- In een for-loop wordt de telling van de basen bijgehouden.
- Een teller wordt geupdated met het aantal reads.
- De phred scores worden berekend per read.
- De langste vijf en de beste vijf reads worden berekend in de best_five() functie. De functie maakt gebruik van een heap-datasctructuur.
- De length en phred distrubutions worden berekend in de functie distrubutions().

Als door elke read van het fastq bestand is gelopen, worden er nog een paar laatste berekeningen gedaan. Zo wordt de gc ratio berekend, de gemiddelde phred en de gemiddelde lengte.
Deze gegevens worden in een JSON file weggeschreven

### filter_trimming.py
Dit script filtert en trimt de reads aan de hand van de instellingen in *config.yaml*. Als input wordt het oude fastq bestand gegeven, de output is het bewerkte fastq bestand.
De functie leest het fastq bestand weer per read in. 
- De sequentie en phred scores worden aan beide kanten getrimd.
- Als de lengte onder de minimale lengte blijkt te zijn, wordt de read overgeslagen en dus ook niet opgeslagen in het nieuwe bestand.
- De gemiddelde phred scores worden berekend. Als deze score onder de minimale phred score valt, wordt deze niet opgeslagen in het nieuwe bestand.

Dit gaat zo door, totdat het fastq bestand helemaal is uitgelezen.

---