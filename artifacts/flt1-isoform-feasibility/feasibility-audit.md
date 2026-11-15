# FLT1 isoform бағыты: feasibility audit

Бұл құжат аннотацияны, processed counts құрылымын және raw reanalysis мүмкіндігін тексереді. PE–FGR айырмашылығы есептелген жоқ.

## Қолдағы annotation

2026-11-14 күні Ensembl GRCh37 REST endpoint-інен алынған `ENSG00000102755` annotation-ында бес protein-coding FLT1 transcript бар. Олар processed GSE114691 кестелерінде де бар. Protein sequence audit төмендегі exact сәйкестіктерді берді:

| Ensembl transcript | UniProt exact match | Атауы |
|---|---|---|
| ENST00000282397 | P17948-1 | Flt1, толық membrane form |
| ENST00000541932 | P17948-3 | sFlt1-14 |
| ENST00000539099 | P17948-4 | isoform 4 |
| ENST00000540678 | P17948-5 | i15 |
| ENST00000543394 | P17948-8 | isoform 8 |

UniProt P17948-2 (`sFlt1`, 687 aa), P17948-6 (`i18`) және P17948-7 (`i21`) осы GRCh37 transcript annotation-ында exact match ретінде жоқ. Сондықтан бес жолдың қосындысын `total soluble FLT1` деп атауға болмайды.

Көздер:

- Ensembl GRCh37 gene lookup: `https://grch37.rest.ensembl.org/lookup/id/ENSG00000102755?expand=1`
- UniProt P17948: `https://www.uniprot.org/uniprotkb/P17948/entry`
- Exact sequence matches: `protein-sequence-matches.json`
- Сақталған Ensembl response: `grch37-flt1-gene.json`

## Transcript-specific exon және junction evidence

GRCh37 annotation-ындағы барлық бес FLT1 transcript бір-біріне қарсы салыстырылды.

| Transcript | Terminal exon, GRCh37 | Басқа FLT1 transcript-терде жоқ contiguous terminal body | Terminal splice junction unique | 50-nt body marker |
|---|---|---:|---|---|
| ENST00000282397, Flt1 | chr13:28874489–28877505 | 2,782 nt | Жоқ; i15 және isoform 8-пен ортақ | Иә |
| ENST00000541932, sFlt1-14 | chr13:28942234–28942800 | 567 nt | Иә | Иә |
| ENST00000539099, isoform 4 | chr13:28973181–28973255 | 75 nt | Иә | Иә, бірақ қысқа |
| ENST00000540678, i15 | chr13:28877281–28877505 | 0 nt | Жоқ; Flt1 және isoform 8-пен ортақ | Жоқ |
| ENST00000543394, isoform 8 | chr13:28877271–28877505 | 0 nt | Жоқ; Flt1 және i15-пен ортақ | Жоқ |

Маңызды айырмашылық: әр terminal exon-ның координатасы бөлек болғанымен, canonical Flt1, i15 және isoform 8 exon-дарының бір splice boundary-і бірдей. Сондықтан толық exon координатасын junction identity ретінде қолдануға болмайды. i15 және isoform 8 terminal body-лары canonical exon ішінде жатыр және terminal junction-дары да ортақ; оларды басқа unique internal exon/junction арқылы ғана ажыратуға болады. Canonical Flt1-дің terminal body sequence-і unique, сондай-ақ unique internal exon-ның екі жағында екі canonical-specific junction бар.

Бұл audit тек FLT1 ішіндегі құрылымдық бірегейлікті көрсетеді. Marker sequence-тердің бүкіл GRCh37 genome/transcriptome бойынша unique mapping-і raw analysis алдында бөлек тексеріледі. Толық кесте: `grch37-flt1-exon-junction-audit.csv`.

Ensembl GRCh37 full cDNA reference-ке жасалған exact 50-mer screen-де primary marker-лер өтті:

- canonical unique internal junction-дардың әрқайсысында 31/31 exact junction-spanning 50-mer unique;
- sFlt1-14 terminal junction-ында 31/31 unique;
- canonical terminal body-де 2,733/2,733 unique;
- sFlt1-14 terminal body-де 493/518 unique.

Бұл transcriptome exact-match screen ғана. Mismatch allowance және genomic alignment behavior таңдалған splice-aware aligner-мен әлі расталуы керек. Script: `audit_50mer_mappability.py`; нәтиже: `marker-50mer-mappability.csv`.

## Processed counts жарамдылығы

GSE114691 processed matrices-тегі бес жол бүтін сандар және бір-бірінен бөлек. Бірақ бұл олардың isoform-specific read assignment-ы сенімді дегенді білдірмейді.

- 79 sample бойынша transcript log-CPM pairwise Spearman correlation-дары 0.919–0.994.
- Әр негізгі transcript жолының бес жол қосындысымен Spearman correlation-ы 0.997–0.998; төмен abundance i15 және isoform 8 үшін 0.933.
- Raw-count үлестері transcript length, shared reads және белгісіз multi-mapping ережесіне тәуелді.
- GEO processing мәтіні Bowtie2 transcriptome alignment және SAM-нан feature abundance алынғанын айтады, бірақ shared reads қалай бөлінгенін, multi-mappers-тің сақталғанын және length/effective-length correction-ды нақтыламайды.

Сондықтан processed count жолдары:

- transcript-тердің reference-де болғанын және signal барын тексеруге жарайды;
- **isoform fraction немесе soluble-to-membrane ratio есептеуге жарамды деп қабылданбайды**;
- PE–FGR contrast үшін қолданылмайды.

Диагностикалық кестелер: `processed-counts-spearman.csv` және `grch37-flt1-exon-junction-audit.csv`.

## Raw data feasibility

GEO GSE114691 raw data-ны SRA `SRP148556` арқылы береді. Ресми run metadata:

- 79 run, бір run бір donor sample;
- барлық кітапхана single-end;
- median read length 50 nt, range 50–58 nt;
- median 41.2 million spots/run;
- барлығы 162.46 Gb sequence және шамамен 62.6 GiB SRA download.

50-nt single-end reads толық isoform reconstruction үшін әлсіз. sFlt1-14 үшін unique terminal body және terminal junction бар. Canonical Flt1 үшін unique terminal body және екі unique internal junction бар, бірақ terminal junction-ның өзі unique емес. i15 және isoform 8 басқа unique internal exon/junction evidence-ке тәуелді және exploratory болып қалады.

Жергілікті дискіде шамамен 12 GiB бос орын бар; толық 79-run reanalysis бұл компьютерде қауіпсіз орындалмайды. Бір run-ды stream/process/delete тәсілімен немесе жеткілікті storage бар external compute-та орындау қажет.

Көздер:

- GEO GSE114691: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691`
- SRA study: `SRP148556`
- Сақталған metadata: `sra-runinfo.csv`, summary: `sra-study-summary.csv`
- ENA FASTQ URL/MD5/size manifest: `ena-fastq-files.tsv` (79 run, шамамен 79.3 GiB compressed FASTQ)

## Шешім

Бағыт **raw-data feasibility stage-ке өтеді**, бірақ processed transcript counts-пен топтық айырмашылық есептелмейді. Annotation және exact cDNA 50-mer gate өтті; splice-aware genome alignment және donor-level raw coverage gate әлі pending.

Primary molecular target ретінде canonical Flt1 (`ENST00000282397`) пен sFlt1-14 (`ENST00000541932`) таңдалады. sFlt1-14 unique terminal junction арқылы, canonical Flt1 unique internal junction-дардың орташа evidence-і арқылы өлшенеді; unique body counts secondary cross-check болады. Қалған үш annotated transcript exploratory feasibility targets ретінде ғана сақталады.

Келесі gate және statistical plan `raw-reanalysis-plan.md` файлында бекітілген. RNA ratio secreted protein мөлшерімен теңестірілмейді және annotation-да жоқ soluble isoform-дар қалпына келтірілді деп айтылмайды.

## Group-blinded raw-read pilot

`SRR7188846` run-ы ENA-дан дискіге сақталмай stream жасалып, алдын ала анықталған exact-unique
50-mer marker-лермен техникалық pilot ретінде тексерілді. Scan кезінде clinical group label қолданылған
жоқ; кейін ресми GEO–SRA manifest бұл run-ның `Control1` екенін растады, сондықтан ол primary PE–FGR
population-ына кірмейді.

- 41,160,188 read толық оқылды; жарамсыз FASTQ record болған жоқ;
- sFlt1-14 unique terminal junction: 101 read;
- canonical-specific internal junction-дар: 4 және 7 read;
- sFlt1-14 және canonical unique terminal body marker-лері: тиісінше 752 және 854 read.

Бұл бір run raw signal барын көрсетеді, бірақ donor coverage gate-ті шешпейді. Атап айтқанда, canonical
junction-дардың біреуі алдын ала бекітілген 5-read шегінен бір read төмен. Сондықтан threshold өзгертілмейді
және бір run нәтижесімен clinical contrast ашылмайды. Exact matching mismatch-терді қабылдамайды әрі
splice-aware genome alignment-ты алмастырмайды. Машина оқитын нәтиже: `raw-pilot-SRR7188846.json`;
scanner: `scan_streamed_fastq_markers.py`.

Official GEO series matrix пен SRA metadata one-to-one join арқылы 79 run-ның manifest-і жасалды:
21 control, 40 PE (`PE only` + `PEIUGR`) және 18 IUGR-only. Файл: `sra-sample-manifest.csv`;
reproduction script: `build_sample_manifest.py`. Бұл mapping coverage/missingness gate үшін қолданылады;
clinical effect gate өтпей ашылмайды.

## Cohort result (2026-11-15, Addendum J1)

**Execution.** All 79 runs were fetched from NCBI SRA, and each run's read count matched the official
`spots` value. FLT1 candidate reads were extracted with `flt1_stream_filter.c` and aligned to the full
GRCh37 genome with HISAT2. Only uniquely aligned reads were counted.

**Aligner mappability check (`hisat2-marker-mappability.json`).** All three primary junction markers
aligned uniquely to the correct position for 31/31 of their 50-mers. Body markers: canonical 0.987,
sFlt1-14 0.868. The body markers are secondary only.

**Pilot regression check.** On `SRR7188846`, the resumable batch version of the aligner reproduced the
single-run counts exactly for every junction marker (113 / 6 / 8). The sFlt1-14 body count differed by
one read (662 vs 661). The likely cause is HISAT2's read-name-seeded tie-breaking, because the batch
version adds a run prefix to read names.

**Group-blinded coverage gate (`cohort-coverage-gate.json`): pass.** Every criterion passed and was
evaluated before any clinical contrast:

- each junction ≥5 reads: 57/58 donors;
- primary junction depth ≥20: 58/58 donors;
- body-ratio vs junction-ratio Spearman rho: 0.879;
- usable donors: 40/40 PE and 17/18 IUGR-only.

Median unique reads per disease donor: sFlt1-14 junction 1,042; canonical junctions 38.5 and 51.5.
About 99.9% of unique FLT1-span reads align to the + strand, while FLT1 lies on the − strand, which is
consistent with a reverse-stranded library.

**Exploratory transcripts.** ENST00000539099, ENST00000540678 and ENST00000543394 had ≥5 junction reads
in 0% of donors. They are reported as not detectable and were not modelled. This detection rule was
committed before the contrast was run.

**Primary contrast (`primary-contrast-results.csv`), run once.**
Model: `junction_ratio ~ PE + log2(FLT1 CPM + 1) + hypoxia_score`, HC3 standard errors, 40 PE vs 18 IUGR-only.

| Analysis | PE coefficient (log2) | 95% CI | p |
|---|---:|---|---:|
| **Primary: junction ratio** | **−0.063** | **−0.524 to 0.398** | **0.79** |
| Secondary: body-read ratio | 0.061 | −0.346 to 0.468 | 0.77 |
| Secondary: PE only vs IUGR only | −0.300 | −0.839 to 0.239 | 0.28 |
| Descriptive: unadjusted | 0.287 | −0.185 to 0.758 | 0.23 |

**Conclusion.** After adjusting for total FLT1 and hypoxia score, this 50-nt single-end dataset shows no
detectable difference between PE and IUGR-only placentas in the sFlt1-14 to canonical Flt1 RNA ratio.
The confidence interval excludes differences larger than roughly ±0.5 log2, i.e. about 1.4-fold. This
does not show that circulating protein is the same in the two groups.

**Descriptive only, not pre-specified.** Both disease groups have higher median junction ratios than
controls: control 3.08, IUGR-only 4.26, PE-only 4.40, PE+IUGR 4.72 (log2). Controls may differ in
gestational age and other factors, so no test was run and this is not a finding.
