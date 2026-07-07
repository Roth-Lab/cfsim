import pandas as pd 

from importlib.resources import files

import cfsim.resources


def run_generate_wig(in_file: str, out_file: str) -> None:

    df_cfdna = (
        pd.read_csv(in_file, sep='\t')
        .assign(
            start=lambda df: df['start'] + 1,     
            end=lambda df: df['end'] + 1
        )
    )
    
    bins_df = bin_hg38(binsize=500000)
    
    df = (
        
        df_cfdna
       
        # get missing bins 
        
        .merge(bins_df, on=['chrom', 'start'], how='outer', suffixes=('_cfdna', '')) 
        
        # fill missing bins
        
        .fillna({'reads': 0})
        
        .astype({'reads': int})
        
        .filter(items=['chrom', 'start', 'end', 'reads'])
        
    )
    
    # convert dataframe to wig
    
    df_to_wig(df=df, genome_build='hg38', output_path=out_file)


def df_to_wig(df: pd.DataFrame, genome_build: str, output_path: str) -> None:
    
    with open(output_path, "w") as f:
        
        for c in df['chrom'].unique():
            
            df_c = (
                
                df
                
                .loc[df['chrom'] == c]
                
                .sort_values("start")
                
            )
            
            start = df_c["start"].iloc[0]
            
            end = df_c['end'].iloc[0]
            
            step = int(end - start)
            
            if start == 0:
                
                start = 1
            
            if genome_build == 'hg19':
                chrom = c.replace('chr')
            else:
                chrom = c
             
            # write header    
            
            s = "fixedStep chrom={chrom} start={start} step={step} span={step}\n".format(
                start=start,
                step=step,
                chrom=chrom
            )
            
            f.write(s)
            
            # write reads
            
            for read in df_c['reads'].values:
                
                f.write("{reads}\n".format(reads=int(read)))



def wig_to_df(wig_file: str) -> pd.DataFrame:
    
    chroms, starts, ends, chrms, values = [], [], [], [], []

    with open(wig_file, "r") as f:
        
        chrom = None
        
        start = None
        
        step = None
        
        span = None
        
        pos = None

        for line in f:
            
            line = line.strip()
            
            if not line:
                
                continue

            if line.startswith("fixedStep"):    # header line 
                
                parts = dict(item.split("=") for item in line.split()[1:])
                
                chrom = parts["chrom"]
                
                start = int(parts["start"])
                
                step = int(parts["step"])
                
                span = int(parts.get("span", step))
                
                pos = start
                
            else:
                
                val = float(line)
                
                chroms.append(chrom)
                
                starts.append(pos)
                
                ends.append(pos + span)
                
                chrms.append(f"{chrom}")
                
                values.append(val)
                
                pos += step
                
    return pd.DataFrame(
        data={
        'chrom': chrms,
        'start': starts,
        'end': ends,
        'values': values
        }
    )



def get_chromosome_sizes() -> pd.DataFrame:
    chrom_sizes_file = files(cfsim.resources).joinpath('hg38.chrom.sizes')
    df = pd.read_csv(chrom_sizes_file, sep='\t', header=None, names=["chrom", "size"])
    chroms = ["chr{}".format(c) for c in range(23)] + ['chrX', 'chrY']
    df = df.loc[df['chrom'].isin(chroms)]
    df['chrom'] = pd.Categorical(df['chrom'], categories=chroms, ordered=True)
    df = df.sort_values(by=['chrom']).reset_index(drop=True)
    return df


def bin_hg38(binsize: int, zero_index: bool = False) -> pd.DataFrame:
    chrom_sizes = get_chromosome_sizes()
    chrom_sizes = chrom_sizes.set_index('chrom')['size'].to_dict()
    rows = []
    for chrom, size in chrom_sizes.items():
        for start in range(0, size, binsize):
            # end = min(start + binsize, size)
            # rows.append((chrom, start, end, end - start))
            rows.append((chrom, start, start + binsize))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    if not zero_index:
        df['start'] += 1
        df['end'] += 1
    return df

