for idcode in 2MQL 5BXA 3RTV 2FNC 2NYP 4Q2Q 2BCQ 5KLC 3BY5 3T6P
do
    # wget https://files.rcsb.org/download/${idcode}.pdb
    sed "s/XXX/${idcode}/" template.dat > parameters.dat
    python3 ~/MMS@FCUL/pypka/pypka parameters.dat
done

