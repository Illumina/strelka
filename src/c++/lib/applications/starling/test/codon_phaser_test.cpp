#include "boost/test/unit_test.hpp"

#include "starling_shared.hh"

#include "blt_util/reference_contig_segment.cpp"
#include "blt_util/seq_util.hh"
#include "blt_common/snp_pos_info.cpp"
#include "starling_common/pos_basecall_buffer.cpp"
#include "starling_common/starling_base_shared.hh"
#include "codon_phaser.cpp"
#include "gvcf_locus_info.cpp"

static void insert_read(const char* read, pos_t position,
        pos_basecall_buffer& bc_buff)
{
    pos_t insert_pos = position;
    while (*read)
    {
        base_call bc(base_to_id(*read), 30, true, 0, 0, false, false, false,
                insert_pos == position, !*(read + 1));
        bc_buff.insert_pos_basecall(insert_pos, true, bc);
        insert_pos++;
        read++;
    }

}

static DIGT::index_t gt_for_base(char base)
{
    static std::map<char, DIGT::index_t> gt =
    {
        { 'A', DIGT::AA },
        { 'C', DIGT::CC },
        { 'G', DIGT::GG },
        { 'T', DIGT::TT }
    };
    return gt[base];
}

BOOST_AUTO_TEST_SUITE( codon_phaser )

BOOST_AUTO_TEST_CASE( simple_3mer )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTAC";
    auto r2 = "ACGTGCTTAC";
    pos_t read_pos = 0;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
        insert_read(r1, read_pos, bc_buff);
    for (int i = 0; i < 10; i++)
        insert_read(r2, read_pos, bc_buff);

    starling_base_options opt;
    opt.phasing_window = 3;

    Codon_phaser phaser(opt, bc_buff);

    std::string phased_variant;

    for (int i = 0; r1[i]; i++)
    {
        site_info si;
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        si.init(read_pos + i, rcs.get_base(read_pos + i), spi, 30);
        si.smod.is_covered = si.smod.is_used_covered = true;
        si.dgt.ref_gt = base_to_id(si.ref);

        if (r1[i] != r2[i])
        {
            si.dgt.is_snp = true;
            si.smod.max_gt = DIGT::AG; // TODO: fix this, but doesn't seem to impact results
        }
        else
        {
            si.smod.max_gt = gt_for_base(r1[i]);
        }
        if (Codon_phaser::is_phasable_site(si) || phaser.is_in_block())
        {
            bool emptyBuffer = phaser.add_site(si);
            if ((!phaser.is_in_block()) || emptyBuffer)
            {
                phased_variant = phaser.buffer()[0].phased_alt;
                phaser.clear();
            }
        }
    }
    BOOST_CHECK_EQUAL("GCT", phased_variant);

}

BOOST_AUTO_TEST_CASE( two_adjacent_3mers )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTACGTACGT";
    auto r2 = "ACGTGCTTGCTTACGT";
    pos_t read_pos = 0;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
        insert_read(r1, read_pos, bc_buff);
    for (int i = 0; i < 10; i++)
        insert_read(r2, read_pos, bc_buff);

    starling_base_options opt;
    opt.phasing_window = 3;

    Codon_phaser phaser(opt, bc_buff);

    std::string phased_variant;

    for (int i = 0; r1[i]; i++)
    {
        site_info si;
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        si.init(read_pos + i, rcs.get_base(read_pos + i), spi, 30);
        si.smod.is_covered = si.smod.is_used_covered = true;
        si.dgt.ref_gt = base_to_id(si.ref);

        if (r1[i] != r2[i])
        {
            si.dgt.is_snp = true;
            si.smod.max_gt = DIGT::AG; // TODO: fix this, but doesn't seem to impact results
        }
        else
        {
            si.smod.max_gt = gt_for_base(r1[i]);
        }
        if (Codon_phaser::is_phasable_site(si) || phaser.is_in_block())
        {
            const bool emptyBuffer = phaser.add_site(si);
            if ((!phaser.is_in_block()) || emptyBuffer)
            {
                phased_variant = phaser.buffer()[0].phased_alt;
                phaser.clear();
            }
        }
    }
    BOOST_CHECK_EQUAL("GCTTGCT", phased_variant);
}

BOOST_AUTO_TEST_CASE( handles_snps_at_start )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTAC";
    auto r2 = "GGGTACGTAC";
    pos_t read_pos = 0;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
        insert_read(r1, read_pos, bc_buff);
    for (int i = 0; i < 10; i++)
        insert_read(r2, read_pos, bc_buff);

    starling_base_options opt;
    opt.phasing_window = 3;

    Codon_phaser phaser(opt, bc_buff);

    std::string phased_variant;

    for (int i = 0; r1[i]; i++)
    {
        site_info si;
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        si.init(read_pos + i, rcs.get_base(read_pos + i), spi, 30);
        si.smod.is_covered = si.smod.is_used_covered = true;
        si.dgt.ref_gt = base_to_id(si.ref);

        if (r1[i] != r2[i])
        {
            si.dgt.is_snp = true;
            si.smod.max_gt = DIGT::AG; // TODO: fix this, but doesn't seem to impact results
        }
        else
        {
            si.smod.max_gt = gt_for_base(r1[i]);
        }
        if (Codon_phaser::is_phasable_site(si) || phaser.is_in_block())
        {
            const bool emptyBuffer = phaser.add_site(si);
            if ((!phaser.is_in_block()) || emptyBuffer)
            {
                phased_variant = phaser.buffer()[0].phased_alt;
                phaser.clear();
            }
        }
    }
    BOOST_CHECK_EQUAL("GG", phased_variant);
}


BOOST_AUTO_TEST_CASE( respects_phasing_window )
{
    reference_contig_segment rcs;
    rcs.seq() = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    pos_basecall_buffer bc_buff(rcs);

    auto r1 = "ACGTACGTACGTACGTACGTACGT";
    auto r2 = "ACGTAAGAACGGAGGTACGTACGT";
    pos_t read_pos = 0;

    // add 2 haplotypes of reads
    for (int i = 0; i < 10; i++)
        insert_read(r1, read_pos, bc_buff);
    for (int i = 0; i < 10; i++)
        insert_read(r2, read_pos, bc_buff);

    starling_base_options opt;
    opt.phasing_window = 5;

    Codon_phaser phaser(opt, bc_buff);

    std::string phased_variant;

    for (int i = 0; r1[i]; i++)
    {
        site_info si;
        const snp_pos_info& spi(bc_buff.get_pos(read_pos + i));
        si.init(read_pos + i, rcs.get_base(read_pos + i), spi, 30);
        si.smod.is_covered = si.smod.is_used_covered = true;
        si.dgt.ref_gt = base_to_id(si.ref);

        if (r1[i] != r2[i])
        {
            si.dgt.is_snp = true;
            si.smod.max_gt = DIGT::AG; // TODO: fix this, but doesn't seem to impact results
        }
        else
        {
            si.smod.max_gt = gt_for_base(r1[i]);
        }
        if (Codon_phaser::is_phasable_site(si) || phaser.is_in_block())
        {
            const bool emptyBuffer = phaser.add_site(si);
            if ((!phaser.is_in_block()) || emptyBuffer)
            {
                phased_variant = phaser.buffer()[0].phased_alt;
                phaser.clear();
            }
        }
    }
    BOOST_CHECK_EQUAL ("AGAACGGAG", phased_variant);
}




BOOST_AUTO_TEST_SUITE_END()

