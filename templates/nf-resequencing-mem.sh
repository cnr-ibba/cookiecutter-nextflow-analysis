# mind to the pipeline version
nextflow run cnr-ibba/nf-resequencing-mem -r v0.6.2 \
    --custom_config_base ${CUSTOM_CONFIG_BASE} \
    -config ${CUSTOM_CONFIG_BASE}/nfcore_custom.config \
    -config ${CUSTOM_CONFIG_BASE}/pipeline/nf-resequencing-mem.config \
    -config conf/custom-{{ cookiecutter.__pipeline }}.config -profile ibba,galileo \
    -resume -params-file conf/params-{{ cookiecutter.__pipeline }}.json
