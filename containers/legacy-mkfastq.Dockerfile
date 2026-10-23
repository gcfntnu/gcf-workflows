# Conversion-only image: never replace the downstream analysis image with this.
# bcl2fastq is statically linked; copy its verified executable without adding
# a second package manager or changing the legacy 10x installation.
ARG BCL2FASTQ_IMAGE=gcfntnu/bcl2fastq@sha256:d63a39ae3ffacb967314a0dfd9b0adb1c52de048aded378e19c32d83588f4ea0
ARG LEGACY_IMAGE=gcfntnu/cellranger@sha256:73fbe65f22fa6274cd1c407a51d7856fbcf9b64990c9c57b427bafb2659951f2
FROM ${BCL2FASTQ_IMAGE} AS bcl2fastq
FROM ${LEGACY_IMAGE}
ARG LEGACY_TOOL=cellranger
COPY --from=bcl2fastq /usr/local/bin/bcl2fastq /usr/local/bin/bcl2fastq
RUN bcl2fastq --version && ${LEGACY_TOOL} --version && ${LEGACY_TOOL} mkfastq --help
