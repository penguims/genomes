# genomes
genomes analysis python scripts

usage: seqman.py [-h]
                 {stat,get,split,cut,translate,transcribe,transform,convert,mutation,chgname}
                 ...

positional arguments:
  {stat,get,split,cut,translate,transcribe,transform,convert,mutation,chgname}
    stat                list sequences stat information
    get                 get a sequence by a name or a list
    split               split a bundle squence to a single sequence file
    cut                 cut a sequence
    translate           translate dna to protein, or in 6 frames
    transcribe          transcribe dna to rna or backward
    transform           reverse and complimant sequence
    convert             convert from one format to another
    mutation            make mutation on sequences
    chgname             change sequence name

optional arguments:
  -h, --help            show this help message and exit

usage: seqman.py stat [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                      [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                      [-al {dna,rna,pro}] [-lo] [-up] [-ds] [-mw] [-3g]
                      [-ck {seguid,gcg,crc32,crc64}] [-ip] [-mp] [-ic]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -ds, --desc           no description when stat
  -mw, --weight         molecular weight
  -3g, --gc123          3 position gc content
  -ck {seguid,gcg,crc32,crc64}, --checksum {seguid,gcg,crc32,crc64}
                        seqeucne checksum
  -ip, --isopoint       isoelectric point
  -mp, --meltpoint      melting point
  -ic, --icc            local component complecity

usage: seqman.py get [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT] [-of OFORMAT]
                     [-fi [FILTER [FILTER ...]]] [-vo] [-al {dna,rna,pro}]
                     [-lo] [-up] [-nm NAME] [-li LIST]
                     [-ft FEATURES [FEATURES ...]]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -nm NAME, --name NAME
                        sequence to get
  -li LIST, --list LIST
                        sequence list
  -ft FEATURES [FEATURES ...], --features FEATURES [FEATURES ...]
                        get sequence by features

usage: seqman.py split [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                       [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                       [-al {dna,rna,pro}] [-lo] [-up] [-di DIR]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -di DIR, --dir DIR    split sequence into dir
  
  usage: seqman.py cut [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT] [-of OFORMAT]
                     [-fi [FILTER [FILTER ...]]] [-vo] [-al {dna,rna,pro}]
                     [-lo] [-up] [-st START] [-en END | -ln LENGTH] [-ra]
                     [-fr FREQUENCY] [-ml MINLEN] [-xl MAXLEN] [-re]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -st START, --start START
                        sequence start
  -en END, --end END    sequence end
  -ln LENGTH, --length LENGTH
                        sequence end
  -ra, --random         random cut a sequence by a length
  -fr FREQUENCY, --frequency FREQUENCY
                        cut times
  -ml MINLEN, --minlen MINLEN
                        mix cut length
  -xl MAXLEN, --maxlen MAXLEN
                        max cut length
  -re, --reverse        cut and concate sequence
  
  sage: seqman.py translate [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                           [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                           [-al {dna,rna,pro}] [-lo] [-up] [-f6]
                           [-ct {Standard,Vertebrate Mitochondrial,Yeast Mitochondrial,Mold Mitochondrial,Protozoan Mitochondrial,Coelenterate Mitochondrial,Mycoplasma; Spiroplasma,Invertebrate Mitochondrial,Ciliate Nuclear,Dasycladacean Nuclear,Hexamita Nuclear,Echinoderm Mitochondrial,Flatworm Mitochondrial,Euplotid Nuclear,Bacterial, Archaeal and Plant Plastid,Alternative Yeast Nuclear,Ascidian Mitochondrial,Alternative Flatworm Mitochondrial,Blepharisma Macronuclear,Chlorophycean Mitochondrial,Trematode Mitochondrial,Scenedesmus obliquus Mitochondrial,Thraustochytrium Mitochondrial,Pterobranchia Mitochondrial,Candidate Division SR1 and Gracilibacteria,Pachysolen tannophilus Nuclear,Karyorelict Nuclear,Condylostoma Nuclear,Mesodinium Nuclear,Peritrich Nuclear,Blastocrithidia Nuclear}]
                           [-lc] [-st START] [-en END] [-fr FRAME]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -f6, --six            translate in 6 frames
  -ct {Standard,Vertebrate Mitochondrial,Yeast Mitochondrial,Mold Mitochondrial,Protozoan Mitochondrial,Coelenterate Mitochondrial,Mycoplasma; Spiroplasma,Invertebrate Mitochondrial,Ciliate Nuclear,Dasycladacean Nuclear,Hexamita Nuclear,Echinoderm Mitochondrial,Flatworm Mitochondrial,Euplotid Nuclear,Bacterial, Archaeal and Plant Plastid,Alternative Yeast Nuclear,Ascidian Mitochondrial,Alternative Flatworm Mitochondrial,Blepharisma Macronuclear,Chlorophycean Mitochondrial,Trematode Mitochondrial,Scenedesmus obliquus Mitochondrial,Thraustochytrium Mitochondrial,Pterobranchia Mitochondrial,Candidate Division SR1 and Gracilibacteria,Pachysolen tannophilus Nuclear,Karyorelict Nuclear,Condylostoma Nuclear,Mesodinium Nuclear,Peritrich Nuclear,Blastocrithidia Nuclear}, --codon {Standard,Vertebrate Mitochondrial,Yeast Mitochondrial,Mold Mitochondrial,Protozoan Mitochondrial,Coelenterate Mitochondrial,Mycoplasma; Spiroplasma,Invertebrate Mitochondrial,Ciliate Nuclear,Dasycladacean Nuclear,Hexamita Nuclear,Echinoderm Mitochondrial,Flatworm Mitochondrial,Euplotid Nuclear,Bacterial, Archaeal and Plant Plastid,Alternative Yeast Nuclear,Ascidian Mitochondrial,Alternative Flatworm Mitochondrial,Blepharisma Macronuclear,Chlorophycean Mitochondrial,Trematode Mitochondrial,Scenedesmus obliquus Mitochondrial,Thraustochytrium Mitochondrial,Pterobranchia Mitochondrial,Candidate Division SR1 and Gracilibacteria,Pachysolen tannophilus Nuclear,Karyorelict Nuclear,Condylostoma Nuclear,Mesodinium Nuclear,Peritrich Nuclear,Blastocrithidia Nuclear}
                        codon table for translation
  -lc, --list           list codon table
  -st START, --start START
                        start translate
  -en END, --end END    start translate
  -fr FRAME, --frame FRAME
                        start translate

usage: seqman.py transcribe [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                            [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                            [-al {dna,rna,pro}] [-lo] [-up]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  
  
  usage: seqman.py transform [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                           [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                           [-al {dna,rna,pro}] [-lo] [-up] [-rv] [-co]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -rv, --reverse        reverse sequence
  -co, --complement     get complememt sequence
  
  usage: seqman.py convert [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                         [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                         [-al {dna,rna,pro}] [-lo] [-up]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence


optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -or ORIGIN, --origin ORIGIN
                        seq id to be changed
  -to TO, --to TO       seq changed to name
  -li LIST, --list LIST
                        sequnce list: from_id to_id
  -kp, --keep           keep origin seq id
  
  sage: seqman.py mutation [-h] [-in INFILE] [-ou OUTPUT] [-if FORMAT]
                          [-of OFORMAT] [-fi [FILTER [FILTER ...]]] [-vo]
                          [-al {dna,rna,pro}] [-lo] [-up]
                          [-tp {snp,del,ins,sv,cnv} [{snp,del,ins,sv,cnv} ...]]
                          [-fr FREQUENCY] [-fc] [-ln LENGTH] [-sp SV_SPAN]
                          [-ra] [-rf] [-ir MIN_RANDOM] [-xr MAX_RANDOM] [-si]
                          [-sy] [-sq SEQ] [-cf CNV_FREQ] [-st START] [-en END]
                          [-rv] [-co] [-rc]

optional arguments:
  -h, --help            show this help message and exit
  -in INFILE, --infile INFILE
                        input file
  -ou OUTPUT, --output OUTPUT
                        output file name
  -if FORMAT, --format FORMAT
                        sequence file format
  -of OFORMAT, --oformat OFORMAT
                        output sequence format
  -fi [FILTER [FILTER ...]], --filter [FILTER [FILTER ...]]
                        sequence filter
  -vo, --verbose        verbose output to stderr
  -al {dna,rna,pro}, --alphabet {dna,rna,pro}
                        sequence alphabet
  -lo, --lower          get lower case sequence
  -up, --upper          get upper case sequence
  -tp {snp,del,ins,sv,cnv} [{snp,del,ins,sv,cnv} ...], --type {snp,del,ins,sv,cnv} [{snp,del,ins,sv,cnv} ...]
                        mutation type
  -fr FREQUENCY, --frequency FREQUENCY
                        mutation frequency
  -fc, --force          force snp to be not same mutation
  -ln LENGTH, --length LENGTH
                        common length to INS, DEL, SV and CNV
  -sp SV_SPAN, --sv_span SV_SPAN
                        SV span length
  -ra, --random         random mutation length
  -rf, --ran_freq       random mutation frequency
  -ir MIN_RANDOM, --min_random MIN_RANDOM
                        min mutation length
  -xr MAX_RANDOM, --max_random MAX_RANDOM
                        max mutation length
  -si, --single         output single mutation sequence
  -sy, --stype          output single mutation type sequence
  -sq SEQ, --seq SEQ    seq to insert, sv or cnv
  -cf CNV_FREQ, --cnv_freq CNV_FREQ
                        cnv repeat times
  -st START, --start START
                        start location for ins, del, sv and cnv
  -en END, --end END    end location for del, sv and cnv
  -rv, --reverse        reverse sequence
  -co, --complement     get complememt sequence
  -rc, --rev_com        get reverse complememt sequence
  
  
