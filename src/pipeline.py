import docker
import os

def seqtk(filepaths,rate):
    list_cmd = []
    for file in filepaths:
        cmd = 'seqtk sample'
        file_ds = file.split('.')[0] + '_ds.fastq'
        cmd += ' ' + file + ' ' + str(rate) +' > ' + file_ds
        list_cmd.append(cmd)
    return list_cmd
def fastqc(filenames):
    cmd = 'fastqc'
    for file in filenames:
        cmd += ' ' + file
    return cmd

def trimmomatic(mode = 'PE',num_threads = 1, p_scale = '-phred33', fwd='',bwd=''):
    cmd = 'trimmomatic'
    if (fwd != '' and bwd != ''):
        cmd += ' ' + mode + ' -threads ' + str(num_threads) + ' ' + p_scale + ' ' + fwd + ' ' + bwd  + ' ' + fwd.split(".")[0] + '_ptrimmed.fastq' + ' ' + fwd.split(".")[0] + '_uptrimmed.fastq' + ' ' + bwd.split(".")[0] + '_trimmed.fastq' + ' ' + bwd.split(".")[0] + '_uptrimmed.fastq' 
    else:
        return 'missing arguments'
    return cmd

def bwa_mem(num_threads, sra_id, platform, ref_genome, fwd, bwd, output):
    cmd = 'bwa mem -t'
    rg_info = '\"@RG\\tID:'+sra_id + '\\tPL:' + platform + '\\tSM:'+ sra_id +'\"'
    cmd += ' '+ str(num_threads) + ' ' + rg_info +  ' ' + ref_genome +  ' ' + fwd +  ' ' + bwd +  ' ' + output

    return cmd


#  ${num_thread} -R ${read_group_info} ${ref} ${reads}/run_forward.fastq.fastq ${reads}/run_backward.fastq ${aligned_reads}/run.paired.bam

def docker_start(dir,docker_dir):
    docker_client = docker.from_env()
    images = docker_client.images.list()
    matching = [img for img in images if "gatk" in img]
    print(matching)
    # matching = ['broadinstitute/gatk:4.4.0.0']
    cmd = 'docker run -v' + ' ' + dir +':' + docker_dir + ' -it' + ' ' + matching[0]
    return cmd

def markDuplicatesSpark(dir, docker_dir):
    docker_start(dir, docker_dir)
    return 1

def main():
    list = ['test/SRR123.fastq', 
            'test/SRR121.fastq',
            'test/SRR120.fastq',
            'test/SRR122.fastq',
            'test/SRR123.fastq',
            'test/SRR124.fastq',
            'test/SRR125.fastq',
            '']
    
    # # fastqc
    # print(fastqc(list))
    # # seqtk
    # for item in seqtk(list, 0.2):
    #     print(item)
    # # trimmomatic
    # for item in list:
    #     print(trimmomatic(fwd=item,bwd=item))
    # print(bwa_mem(1,'SRR2020636','ILLUMINA','test/hg38.fa',list[0], list[1], 'paired.bam'))
    # print(docker_start('data','/gatk/data'))
    markDuplicatesSpark('.','/gatk/demo/')

if __name__ == "__main__":
    main()