#ifndef DSPCONFIG_H_
#define DSPCONFIG_H_

/** The version should be synced with the one in Github.
 * So it is based on symantic versioning.
 */
#define DSP_VERSION_MAJOR 1
#define DSP_VERSION_MINOR 4
#define DSP_VERSION_PATCH 1

#include <stdio.h>

inline void show_copyright()
{
    char msg[1024];
    sprintf(msg, "\n=================================================================================\n"
                 "  DSP: Parallel decomposition methods for structured programming\n"
                 "  - Version %d.%d.%d\n"
                 "  - See https://github.com/Argonne-National-Laboratory/DSP\n\n"
                 "  Under the terms of Contract No. DE-AC02-06CH11357 with UChicago Argonne, LLC,\n"
                 "  the U.S. Government retains certain rights in this software.\n\n"
                 "=================================================================================\n",
            DSP_VERSION_MAJOR, DSP_VERSION_MINOR, DSP_VERSION_PATCH);
    printf("%s", msg);
}

#endif
