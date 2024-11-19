import os

import dt4dds

# set up config
dt4dds.default_logging()
dt4dds.config.show_progressbars = True


# define primer sequences for PCR定义PCR的引物序列,默认为后期添加引物，无需自己设计
primers_0 = ["ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT"]
primers_2 = ["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"]
#input_file = '/mnt/00.zhibo/random/row_min.fasta'
#output_file= "/mnt/00.zhibo/random/"


def chuanxing(input_row_fasta,output_file,cycle=1):

    # assign efficiency properties for amplification, here a normal distribution with std of 0.0051为扩增分配效率属性，这里是 std 为 0.0051 的正态分布
    dt4dds.properties.set_property_settings(
        dt4dds.settings.defaults.SequenceProperties(
            efficiency_distribution='uniform',
            efficiency_params={'low':0.8 , 'high': 1.2},
        )
    )

    # read the sequences from the provided example file读取提供的例子
    seq_list = dt4dds.tools.txt_to_seqlist(input_row_fasta)  # ("design_sequences.txt")
    n_seqs = len(seq_list)

    # 使用电化学合成的默认值设置
    synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Photolithographic()
    # settings can be customized further when passing to the process instance
    array_synthesis = dt4dds.processes.ArraySynthesis(synthesis_settings(
        oligo_distribution_type='lognormal',
        oligo_distribution_params={'mean': 0, 'sigma': 0.30},
    ))
    array_synthesis.process(seq_list)
    # print("array",array_synthesis.process(seq_list))
    # print('nseq',n_seqs)
    # sample with mean coverage of 200平均覆盖率为 200 的样本,可以按重量或体积对池进行采样
    pool = array_synthesis.sample_by_counts(500 * n_seqs)
    # print("pool1========\n",pool)
    pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
    # pool.volume = 1

    #
    # Aging for one half-live使用半衰期
    #
    aging_settings = dt4dds.settings.defaults.Aging()
    aging = dt4dds.processes.Aging(aging_settings(), n_halflives=1)
    pool = aging.process(pool)
    pool.volume = 1

    for i in range(cycle):
        if i > 0:
            name = 'chuan'
        else:
            name = ''
        output_filename = os.path.join(output_file,f"ROW_{name}R{i}")
        # 使用高保真聚合酶进行 20 个循环的 PCR，平均效率为 95%
        pcr_settings = dt4dds.settings.defaults.PCR_HiFi()
        pcr = dt4dds.processes.PCR(pcr_settings(
            primers=primers_0,
            template_volume=1,
            volume=20,
            efficiency_mean=0.95,
            n_cycles=20,
        ))
        pool = pcr.process(pool)
        # Sequencing-by-synthesis with paired reads and sequencing coverage of 25
        # 序列合成配对读取和测序覆盖率为 25
        synthesis_settings = dt4dds.settings.defaults.iSeq100()
        sbs_sequencing = dt4dds.processes.SBSSequencing(synthesis_settings(
            output_directory=output_filename,#Create a new folder
            n_reads=int(200 * n_seqs),
            read_length=150,
            read_mode='single-end',
        ))

        #ues this to chuan
        pool = sbs_sequencing.process(pool)
        pool.volume = 1
        print("end----------------------")

def bingxing(input_row_fasta,output_file,cycle=1):

    # 使用高保真聚合酶进行 30 个循环的 PCR，平均效率为 95%
    for i in range(cycle):
        # assign efficiency properties for amplification, here a normal distribution with std of 0.0051为扩增分配效率属性，这里是 std 为 0.0051 的正态分布
        dt4dds.properties.set_property_settings(
            dt4dds.settings.defaults.SequenceProperties(
                efficiency_distribution='normal',
                efficiency_params={'loc': 1.0, 'scale': 0.0051},
            )
        )

        # read the sequences from the provided example file读取提供的例子
        seq_list = dt4dds.tools.txt_to_seqlist(input_row_fasta)  # ("design_sequences.txt")
        n_seqs = len(seq_list)

        # 使用电化学合成的默认值设置
        synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Twist()
        # settings can be customized further when passing to the process instance
        array_synthesis = dt4dds.processes.ArraySynthesis(synthesis_settings(
            oligo_distribution_type='lognormal',
            oligo_distribution_params={'mean': 0, 'sigma': 0.30},
        ))
        array_synthesis.process(seq_list)
        # print("array",array_synthesis.process(seq_list))
        # print('nseq',n_seqs)
        # sample with mean coverage of 200平均覆盖率为 200 的样本,可以按重量或体积对池进行采样
        pool = array_synthesis.sample_by_counts(500 * n_seqs)
        # print("pool1========\n",pool)
        pool = dt4dds.generators.attach_primers_to_pool(pool, *primers_0)
        # pool.volume = 1

        #
        # Aging for one half-live使用半衰期
        #
        aging_settings = dt4dds.settings.defaults.Aging()
        aging = dt4dds.processes.Aging(aging_settings(), n_halflives=1)
        pool = aging.process(pool)
        pool.volume = 1

        output_filename = os.path.join(output_file, f"ROW_bingR{i}")
        pcr_settings = dt4dds.settings.defaults.PCR_HiFi()
        pcr = dt4dds.processes.PCR(pcr_settings(
            primers=primers_0,
            template_volume=1,
            volume=20,
            efficiency_mean=0.95,
            n_cycles=20,
        ))
        pool = pcr.process(pool)
        # Sequencing-by-synthesis with paired reads and sequencing coverage of 25
        # 序列合成配对读取和测序覆盖率为 25
        synthesis_settings = dt4dds.settings.defaults.iSeq100()
        sbs_sequencing = dt4dds.processes.SBSSequencing(synthesis_settings(
            output_directory=output_filename,#Create a new folder
            n_reads=int(200 * n_seqs),
            read_length=150,
            read_mode='single-end',
        ))
        # single
        sbs_sequencing.process(pool)
        print("end----------------------")


