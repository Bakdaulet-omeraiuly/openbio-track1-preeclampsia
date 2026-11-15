# FLT1 isoform raw RNA-seq reanalysis жоспары

Бұл жоспар clinical-group contrast қаралмай тұрып жазылды. Feasibility gate өтпесе, PE–FGR тесті жүргізілмейді.

## Зерттеу сұрағы

Жалпы FLT1 activity және placental hypoxia score ескерілгенде, PE бар жүктіліктерде sFlt1-14-тің canonical membrane Flt1-ге қатысты RNA evidence-і FGR-only жүктіліктерден өзгеше ме?

## Нысаналар

- **Primary numerator:** sFlt1-14, ENST00000541932, unique terminal junction `chr13:28942801-28959021` және unique terminal exon body.
- **Primary denominator:** canonical Flt1, ENST00000282397, екі unique internal junction-ның орташа count-ы (`chr13:28931823-28959021` және `chr13:28919689-28931690`). Canonical terminal junction i15/isoform 8-пен ортақ болғандықтан қолданылмайды; unique terminal body secondary cross-check болады.
- **Exploratory:** ENST00000539099, ENST00000540678 және ENST00000543394. i15/isoform 8 terminal body-лары canonical exon ішінде жатқандықтан, junction немесе басқа unique exon evidence қажет.
- `Total soluble FLT1 fraction` есептелмейді: GRCh37 annotation UniProt sFlt1 P17948-2, i18 және i21 түрлерін толық қамтымайды.

## Reference және alignment

1. GRCh37 genome FASTA және audit-пен бір release-ке сәйкес Ensembl GTF қолданылады.
2. Synthetic 50-mer және junction-spanning sequences толық genome/transcriptome-ға салыстырылады. Primary marker-лердің кемінде бір uniquely mappable body немесе junction evidence-і болуы міндетті.
3. Raw single-end reads бір run-нан бір уақытта stream/process етіледі; аралық FASTQ өңдеуден кейін жойылады.
4. Splice-aware genomic alignment `sjdbOverhang=49` параметрімен орындалады. Strandedness metadata-да анық емес; gene-body және splice evidence арқылы тексеріліп, analysis log-қа жазылады.
5. Primary evidence тек uniquely mapped reads-тен алынады:
   - exact terminal junction-spanning reads;
   - басқа FLT1 transcript exon-дарымен қабыспайтын terminal body reads.
6. Transcript-level probabilistic quantification, егер орындалса, secondary sensitivity болады. 50-nt single-end fragment-length uncertainty және bootstrap uncertainty міндетті түрде сақталады.

## Feasibility gate, group labels қолданылмайды

Primary clinical model тек төмендегі шарттардың бәрі орындалса іске қосылады:

1. Canonical және sFlt1-14 target-тарының әрқайсысында толық reference бойынша unique mapping marker бар.
2. 58 disease donor-дың кемінде 70%-ында sFlt1-14 junction және canonical junction-дардың әрқайсысы бойынша кемінде 5 unique read бар.
3. Donor-лардың кемінде 70%-ында canonical + sFlt1-14 бойынша жиынтық primary junction depth кемінде 20 read.
4. Length-normalized unique-body ratio мен junction ratio donor деңгейінде Spearman rho ≥ 0.70 немесе олардың айырмашылығы техникалық себеппен түсіндіріліп, бір өлшем primary ретінде алдын ала бекітіледі.
5. Бір clinical group-та usable sample үлесі 60%-дан төмен болмейді. Бұл check group effect-ті қарамай, тек missingness/coverage-ті тексереді.

Gate өтпесе, нәтиже `available 50-nt single-end data do not support a reliable isoform-composition test` деп жабылады. Processed count rows substitute ретінде қолданылмайды.

Status: exact 50-mer screen against the full Ensembl GRCh37 cDNA reference passed for all predefined primary markers. The same markers still require confirmation with the selected splice-aware genome aligner before donor coverage is evaluated.

Technical pilot status: one streamed run (`SRR7188846`, scanned without its group label and later mapped to
`Control1`) contained 101 exact reads for the
sFlt1-14 junction and 4/7 exact reads for the two canonical-specific junctions. This establishes detectable
raw signal but does not pass or fail the cohort gate. The pre-specified thresholds above are unchanged.

## Primary outcome және contrast

Gate өткен жағдайда:

- Donor-level primary outcome: `log2((sFlt1-14 unique terminal-junction reads + 0.5) / (mean of the two canonical unique internal-junction counts + 0.5))`.
- Primary population: 40 PE donors (`PEONLY` + `PEIUGR`) versus 18 `IUGRONLY` donors.
- Primary model: `isoform_log_ratio ~ PE + log2(total_FLT1_CPM + 1) + hypoxia_score` with HC3 robust standard errors.
- Primary test: two-sided PE coefficient, alpha 0.05.
- Controls are descriptive context only and do not enter the primary PE–FGR contrast.
- Secondary sensitivity: body-read ratio; PE-only versus IUGR-only; probabilistic transcript abundance if quantification uncertainty is acceptable.
- Other three transcripts are exploratory; their p-values receive Benjamini–Hochberg correction as one family.

## Interpretation

- Positive result: evidence for a difference in **RNA isoform composition** in this dataset, not proof of different circulating sFlt-1 protein or secretion.
- Null result: no detectable composition difference at the resolution and coverage of 50-nt single-end RNA-seq; it does not prove isoforms or protein effects are identical.
- No result is called novel until a targeted literature review confirms that the exact PE-versus-FGR isoform contrast has not already been reported.

## Compute requirement

SRA archive size is about 62.6 GiB; the official ENA manifest totals about 79.3 GiB of compressed FASTQ,
and uncompressed reads are much larger. The current machine has insufficient free disk for a full retained
download. Use external storage/compute or a one-run-at-a-time streaming workflow with the saved ENA URLs,
MD5 checksums, run-level logs and resumable outputs.

The pure-Python exact-marker scanner took roughly 14 minutes for a 41.2-million-read pilot and is retained
only for audit/pilot reproduction. Cohort-scale production should use a compiled splice-aware aligner and
compiled junction counter on external storage/compute.
