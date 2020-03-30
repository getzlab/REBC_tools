FROM broadgdac/matlab:2016a
## Install Xvfb
RUN set -ex \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        xvfb \
        xfonts-100dpi \
        xfonts-75dpi \
        xfonts-cyrillic \
        x-ttcidfont-conf \
    && rm -rf /var/lib/apt/lists/*

COPY src /src
WORKDIR src
RUN chmod +x /src/fc_MAF_AC_PP_CCF_fit_v3 /src/run_fc_MAF_AC_PP_CCF_fit_v3.sh; sync;
