#paste in order of phylogeny, unsure if need but was nervous
paste ortho.cap.hypot.tpm.tsv ortho.ovi.hypot.tpm.tsv ortho.bos.hypot.tpm.tsv ortho.sus.hypot.tpm.tsv ortho.sor.hypot.tpm.tsv ortho.cav.hypot.tpm.tsv ortho.het.hypot.tpm.tsv ortho.fuk.hypot.tpm.tsv ortho.per.hypot.tpm.tsv ortho.mic.hypot.tpm.tsv ortho.mus.hypot.tpm.tsv ortho.rat.hypot.tpm.tsv ortho.ict.hypot.tpm.tsv ortho.mac.hypot.tpm.tsv ortho.hom.hypot.tpm.tsv ortho.pan.hypot.tpm.tsv > temp2
### sorry for this long awk script, but EVE is tempermental with format.
### print $11 first to get the ids first (everything in terms of bos)

awk '{print $11 FS $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $14 FS $15 FS $16 FS $17 FS $18 FS $19 FS $20 FS $21 FS $22 FS $23 FS $24 FS $25 FS $26 FS $27 FS $28 FS $29 FS $30 FS $31 FS $32 FS $33 FS $34 FS $35 FS $36 FS $37 FS $38 FS $39 FS $40 FS $41 FS $42 FS $43 FS $44 FS $45 FS $46 FS $47 FS $48 FS $49 FS $50 FS $51 FS $52 FS $53 FS $54 FS $55 FS $56 FS $57 FS $58 FS $59 FS $60 FS $61 FS $62 FS $63 FS $64 FS $65 FS $66 FS $67 FS $68 FS $69 FS $70 FS $71 FS $72 FS $73 FS $74 FS $75 FS $76 FS $77 FS $78 FS $79 FS $80 FS $81 FS $82 FS $83 FS $84 FS $85 FS $86 FS $87 FS $88 FS $89 FS $90 FS $91 FS $92 FS $93 FS $94 FS $95 FS $96 FS $97 FS $98}' temp2 > adult_nopap_hypothALV