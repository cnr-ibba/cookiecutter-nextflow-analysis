# calling the pipeline
nextflow run nf-core/sarek -r 3.8.1 \
    --custom_config_base ${CUSTOM_CONFIG_BASE} \
    -config ${CUSTOM_CONFIG_BASE}/nfcore_custom.config \
    -config ${CUSTOM_CONFIG_BASE}/pipeline/sarek.config \
    -config conf/custom-{{ cookiecutter.__pipeline }}.config -profile ibba,galileo \
    -resume -params-file conf/params-{{ cookiecutter.__pipeline }}.json
