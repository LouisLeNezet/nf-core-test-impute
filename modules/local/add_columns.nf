process ADD_COLUMNS {
    label 'process_single'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'NR==2{
        print \$0, "best_gt_rsquared", "imputed_ds_rsquared", "ID", "Region", "Depth", "Panel"
    } NR>2 && NR<=10 {
        \$(NF+1)="${meta.id}"
        \$(NF+1)="${meta.region}"
        \$(NF+1)="${meta.depth}"
        \$(NF+1)="${meta.panel}"
        print
    }' $input > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -1 | grep -o -E '([0-9]+.){1,2}[0-9]')
    END_VERSIONS
    """
}
