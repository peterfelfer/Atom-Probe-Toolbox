#!/bin/bash
# Fix 143 broken symlinks in exampleFiles
# These files only exist locally (not in Dropbox), so we link to the backup.
# Run this from the Atom-Probe-Toolbox directory on your Mac.

TOOLBOX_DIR="$(cd "$(dirname "$0")" && pwd)"
echo "Fixing broken symlinks in $TOOLBOX_DIR/exampleFiles..."

fixed=0
failed=0

if [ -L "$TOOLBOX_DIR/exampleFiles/R6025_267386/R6025_267386_NMC811.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R6025_267386/R6025_267386_NMC811.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R6025_267386_NMC811.pos" "$TOOLBOX_DIR/exampleFiles/R6025_267386/R6025_267386_NMC811.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/77c8293e-1b1c-41ca-96c9-e18ff266a000.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/77c8293e-1b1c-41ca-96c9-e18ff266a000.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/77c8293e-1b1c-41ca-96c9-e18ff266a000.POS" "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/77c8293e-1b1c-41ca-96c9-e18ff266a000.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/77c8293e-1b1c-41ca-96c9-e18ff266a000/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02_xrng.fig" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02_xrng.fig" && \
  ln -s "../exampleFiles_backup/castrip grain boundary/R04_13509-v02_xrng.fig" "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02_xrng.fig" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.pos" && \
  ln -s "../exampleFiles_backup/castrip grain boundary/R04_13509-v02.pos" "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.xrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.xrng" && \
  ln -s "../exampleFiles_backup/castrip grain boundary/R04_13509-v02.xrng" "$TOOLBOX_DIR/exampleFiles/R04_13509/R04_13509-v02.xrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/8fe936a9-c8bf-417a-8006-bff810375bc7.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/8fe936a9-c8bf-417a-8006-bff810375bc7.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/8fe936a9-c8bf-417a-8006-bff810375bc7.POS" "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/8fe936a9-c8bf-417a-8006-bff810375bc7.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/8fe936a9-c8bf-417a-8006-bff810375bc7/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471.POS" "$TOOLBOX_DIR/exampleFiles/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471/ee04a8e4-2f1a-472d-89a0-37a0a3fd6471.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01/Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R2001_00016/R2001_00016 NiAlFeCo.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R2001_00016/R2001_00016 NiAlFeCo.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R2001_00016 NiAlFeCo.pos" "$TOOLBOX_DIR/exampleFiles/R2001_00016/R2001_00016 NiAlFeCo.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/39f09588-44a8-4dc8-862b-1db9790ee30f.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/39f09588-44a8-4dc8-862b-1db9790ee30f.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/39f09588-44a8-4dc8-862b-1db9790ee30f.POS" "$TOOLBOX_DIR/exampleFiles/39f09588-44a8-4dc8-862b-1db9790ee30f/39f09588-44a8-4dc8-862b-1db9790ee30f.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01.epos" && \
  ln -s "../exampleFiles_backup/Si/R04_19434-v01.epos" "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01_rng.fig" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01_rng.fig" && \
  ln -s "../exampleFiles_backup/Si/R04_19434-v01_rng.fig" "$TOOLBOX_DIR/exampleFiles/R04_19434/R04_19434-v01_rng.fig" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03/Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_50673/R18_50673-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_50673/R18_50673-v01.pos" && \
  ln -s "../exampleFiles_backup/unranged/2013ESI_nc_A220/R18_50673-v01.pos" "$TOOLBOX_DIR/exampleFiles/R18_50673/R18_50673-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Pristine FeO R76_45587/Pristine FeO R76_45587.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Pristine FeO R76_45587/Pristine FeO R76_45587.EPOS" && \
  ln -s "../exampleFiles_backup/external_repos/5838237/Pristine FeO R76_45587.EPOS" "$TOOLBOX_DIR/exampleFiles/Pristine FeO R76_45587/Pristine FeO R76_45587.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/bb158535-0204-4619-b532-e35c422c96b5.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/bb158535-0204-4619-b532-e35c422c96b5.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/bb158535-0204-4619-b532-e35c422c96b5.POS" "$TOOLBOX_DIR/exampleFiles/bb158535-0204-4619-b532-e35c422c96b5/bb158535-0204-4619-b532-e35c422c96b5.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc.POS" "$TOOLBOX_DIR/exampleFiles/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc/653b9b5f-6f4c-4773-8e80-e0e9f3fe5dcc.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/485df44e-33bd-4c40-89fb-29aa3faa7785.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/485df44e-33bd-4c40-89fb-29aa3faa7785.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/485df44e-33bd-4c40-89fb-29aa3faa7785.POS" "$TOOLBOX_DIR/exampleFiles/485df44e-33bd-4c40-89fb-29aa3faa7785/485df44e-33bd-4c40-89fb-29aa3faa7785.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m14 Young OES850 --19M R31_20642-v01.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m14 Young OES850 --19M R31_20642-v01.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m14 Young OES850 --19M R31_20642-v01/Enamel-2 m14 Young OES850 --19M R31_20642-v01.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R6012_2640257000/R6012_2640257000 Al.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R6012_2640257000/R6012_2640257000 Al.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R6012_2640257000 Al.pos" "$TOOLBOX_DIR/exampleFiles/R6012_2640257000/R6012_2640257000 Al.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0.POS" "$TOOLBOX_DIR/exampleFiles/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0/6c6c6ee5-4080-46bf-8cfb-d33c0863fdd0.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/70_50_50/AB.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/70_50_50/AB.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/APM.LEAP.Datasets.4 2/AB.RRNG" "$TOOLBOX_DIR/exampleFiles/70_50_50/AB.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/70_50_50/70_50_50.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/70_50_50/70_50_50.POS" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/APM.LEAP.Datasets.4/70_50_50.POS" "$TOOLBOX_DIR/exampleFiles/70_50_50/70_50_50.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v06.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v06.EPOS" && \
  ln -s "../exampleFiles_backup/R56_08424_FHT_NiAl/R56_08424-v06.EPOS" "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v06.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v04.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v04.pos" && \
  ln -s "../exampleFiles_backup/R56_08424_FHT_NiAl/R56_08424-v04.pos" "$TOOLBOX_DIR/exampleFiles/R56_08424/R56_08424-v04.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/1min FeO R76_50315/1min FeO R76_50315.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/1min FeO R76_50315/1min FeO R76_50315.EPOS" && \
  ln -s "../exampleFiles_backup/external_repos/5838237/1min FeO R76_50315.EPOS" "$TOOLBOX_DIR/exampleFiles/1min FeO R76_50315/1min FeO R76_50315.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.xrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.xrng" && \
  ln -s "../exampleFiles_backup/Alex_data/R18_52010-v02.xrng" "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.xrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.pos" && \
  ln -s "../exampleFiles_backup/Alex_data/R18_52010-v02.pos" "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02_xrng.fig" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02_xrng.fig" && \
  ln -s "../exampleFiles_backup/Alex_data/R18_52010-v02_xrng.fig" "$TOOLBOX_DIR/exampleFiles/R18_52010/R18_52010-v02_xrng.fig" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_51491/2013-11-04_51491-V14.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_51491/2013-11-04_51491-V14.RRNG" && \
  ln -s "../exampleFiles_backup/2014Pt triple_line_Scherrer/2013-11-04_51491-V14.RRNG" "$TOOLBOX_DIR/exampleFiles/R18_51491/2013-11-04_51491-V14.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_51491/R18_51491-v14.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_51491/R18_51491-v14.pos" && \
  ln -s "../exampleFiles_backup/2014Pt triple_line_Scherrer/R18_51491-v14.pos" "$TOOLBOX_DIR/exampleFiles/R18_51491/R18_51491-v14.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db.POS" "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/3a5ddf5d-ddf2-4c8c-89fe-e1dd115683db/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m17 Old OES850 --17M R31_20647-v01.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m17 Old OES850 --17M R31_20647-v01.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m17 Old OES850 --17M R31_20647-v01/Enamel-2 m17 Old OES850 --17M R31_20647-v01.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R13_40310/R13_40310Zr Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R13_40310/R13_40310Zr Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R13_40310Zr Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R13_40310/R13_40310Zr Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/c55d7256-3600-4e0a-91d2-f09fd4648ef9.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/c55d7256-3600-4e0a-91d2-f09fd4648ef9.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/c55d7256-3600-4e0a-91d2-f09fd4648ef9.POS" "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/c55d7256-3600-4e0a-91d2-f09fd4648ef9.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/c55d7256-3600-4e0a-91d2-f09fd4648ef9/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R16_50678/R16_50678 NiPdSi Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R16_50678/R16_50678 NiPdSi Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R16_50678 NiPdSi Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R16_50678/R16_50678 NiPdSi Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/d5692e46-52a1-4c62-a0b2-01fe39c2a496.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/d5692e46-52a1-4c62-a0b2-01fe39c2a496.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/d5692e46-52a1-4c62-a0b2-01fe39c2a496.POS" "$TOOLBOX_DIR/exampleFiles/d5692e46-52a1-4c62-a0b2-01fe39c2a496/d5692e46-52a1-4c62-a0b2-01fe39c2a496.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5076_69138/R5076_69138-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5076_69138/R5076_69138-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/10677562/R5076_69138-v01.epos" "$TOOLBOX_DIR/exampleFiles/R5076_69138/R5076_69138-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R76_23219/R76_23219-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R76_23219/R76_23219-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/2669482/R76_23219-v01.epos" "$TOOLBOX_DIR/exampleFiles/R76_23219/R76_23219-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5100_235274/R5100_235274UW Ranging - Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5100_235274/R5100_235274UW Ranging - Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R5100_235274UW Ranging - Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R5100_235274/R5100_235274UW Ranging - Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R17_89420/R17_89420FeSuper.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R17_89420/R17_89420FeSuper.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R17_89420FeSuper.pos" "$TOOLBOX_DIR/exampleFiles/R17_89420/R17_89420FeSuper.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_53222/R18_53222_W_18K-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_53222/R18_53222_W_18K-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/R18_53222_W_18K-v01.epos" "$TOOLBOX_DIR/exampleFiles/R18_53222/R18_53222_W_18K-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19410/R04_19410-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19410/R04_19410-v01.epos" && \
  ln -s "../exampleFiles_backup/unranged/Al/R04_19410-v01.epos" "$TOOLBOX_DIR/exampleFiles/R04_19410/R04_19410-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01/Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R15_72243/R15_72243 Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R15_72243/R15_72243 Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R15_72243 Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R15_72243/R15_72243 Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R6006_254276/R6006_254276 Cu VpL.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R6006_254276/R6006_254276 Cu VpL.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R6006_254276 Cu VpL.pos" "$TOOLBOX_DIR/exampleFiles/R6006_254276/R6006_254276 Cu VpL.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.pos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_portland_wang/R31_06365-v02.pos" "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/APM.LEAP.Datasets.1 2/R31_06365-v02.rrng" "$TOOLBOX_DIR/exampleFiles/R31_06365/R31_06365-v02.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19421/R04_19421-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19421/R04_19421-v01.epos" && \
  ln -s "../exampleFiles_backup/unranged/W/R04_19421-v01.epos" "$TOOLBOX_DIR/exampleFiles/R04_19421/R04_19421-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.pos" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m19 old OES850 --11M R31_20650-v03.pos" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/11111719/Enamel-2 m19 old OES850 --11M R31_20650-v03.rrng" "$TOOLBOX_DIR/exampleFiles/Enamel-2 m19 old OES850 --11M R31_20650-v03/Enamel-2 m19 old OES850 --11M R31_20650-v03.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/ee98f685-6f6f-484d-b8ea-9980930f174d.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/ee98f685-6f6f-484d-b8ea-9980930f174d.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/ee98f685-6f6f-484d-b8ea-9980930f174d.POS" "$TOOLBOX_DIR/exampleFiles/ee98f685-6f6f-484d-b8ea-9980930f174d/ee98f685-6f6f-484d-b8ea-9980930f174d.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R40_110290/R40_110290Au.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R40_110290/R40_110290Au.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R40_110290Au.pos" "$TOOLBOX_DIR/exampleFiles/R40_110290/R40_110290Au.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5038_00333/rng_5pj.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5038_00333/rng_5pj.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_gbco/rng_5pj.rrng" "$TOOLBOX_DIR/exampleFiles/R5038_00333/rng_5pj.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5038_00333/R5038_00333-v02.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5038_00333/R5038_00333-v02.epos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_gbco/R5038_00333-v02.epos" "$TOOLBOX_DIR/exampleFiles/R5038_00333/R5038_00333-v02.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R76_20231/R76_20231-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R76_20231/R76_20231-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/2669484/R76_20231-v01.epos" "$TOOLBOX_DIR/exampleFiles/R76_20231/R76_20231-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/70a59eff-003c-4337-832a-604c260dc623.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/70a59eff-003c-4337-832a-604c260dc623.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/70a59eff-003c-4337-832a-604c260dc623.POS" "$TOOLBOX_DIR/exampleFiles/70a59eff-003c-4337-832a-604c260dc623/70a59eff-003c-4337-832a-604c260dc623.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01_heat treated.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01_heat treated.pos" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/R02_10806-v01_heat treated.pos" "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01_heat treated.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01.pos" && \
  ln -s "../exampleFiles_backup/unranged/ReneN5-HT/R02_10806-v01.pos" "$TOOLBOX_DIR/exampleFiles/R02_10806/R02_10806-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_10806/Rene N5.RNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_10806/Rene N5.RNG" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/Rene N5.RNG" "$TOOLBOX_DIR/exampleFiles/R02_10806/Rene N5.RNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680-v02.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680-v02.pos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/aut_leoben_leitner/R21_08680-v02.pos" "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680-v02.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680.rrng" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/aut_leoben_leitner/R21_08680.rrng" "$TOOLBOX_DIR/exampleFiles/R21_08680/R21_08680.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R21_07055/7055_range file_analysis.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R21_07055/7055_range file_analysis.rrng" && \
  ln -s "../exampleFiles_backup/2015Visit_Kathi/7055_range file_analysis.rrng" "$TOOLBOX_DIR/exampleFiles/R21_07055/7055_range file_analysis.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R21_07055/R21_07055-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R21_07055/R21_07055-v01.pos" && \
  ln -s "../exampleFiles_backup/2015Visit_Kathi/R21_07055-v01.pos" "$TOOLBOX_DIR/exampleFiles/R21_07055/R21_07055-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5076_69145/R5076_69145-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5076_69145/R5076_69145-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/10677562/R5076_69145-v01.epos" "$TOOLBOX_DIR/exampleFiles/R5076_69145/R5076_69145-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/5acb03b0-04db-4a1d-93a1-356cbd913e69.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/5acb03b0-04db-4a1d-93a1-356cbd913e69.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/5acb03b0-04db-4a1d-93a1-356cbd913e69.POS" "$TOOLBOX_DIR/exampleFiles/5acb03b0-04db-4a1d-93a1-356cbd913e69/5acb03b0-04db-4a1d-93a1-356cbd913e69.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/b4b7c711-c302-4ee5-b411-e2a3bf6aad97.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/b4b7c711-c302-4ee5-b411-e2a3bf6aad97.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/b4b7c711-c302-4ee5-b411-e2a3bf6aad97.POS" "$TOOLBOX_DIR/exampleFiles/b4b7c711-c302-4ee5-b411-e2a3bf6aad97/b4b7c711-c302-4ee5-b411-e2a3bf6aad97.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/dd15537d-6019-4a80-8abe-4b388f3b77b6.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/dd15537d-6019-4a80-8abe-4b388f3b77b6.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/dd15537d-6019-4a80-8abe-4b388f3b77b6.POS" "$TOOLBOX_DIR/exampleFiles/dd15537d-6019-4a80-8abe-4b388f3b77b6/dd15537d-6019-4a80-8abe-4b388f3b77b6.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_13447/R02_13447-v01_Rafted.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_13447/R02_13447-v01_Rafted.pos" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/R02_13447-v01_Rafted.pos" "$TOOLBOX_DIR/exampleFiles/R02_13447/R02_13447-v01_Rafted.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_13447/Rene N5.RNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_13447/Rene N5.RNG" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/Rene N5.RNG" "$TOOLBOX_DIR/exampleFiles/R02_13447/Rene N5.RNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_52189/2013-12-16_R18_52189-v01.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_52189/2013-12-16_R18_52189-v01.rrng" && \
  ln -s "../exampleFiles_backup/140618KIT_anode/2013-12-16_R18_52189-v01.rrng" "$TOOLBOX_DIR/exampleFiles/R18_52189/2013-12-16_R18_52189-v01.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_52189/R18_52189-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_52189/R18_52189-v01.epos" && \
  ln -s "../exampleFiles_backup/140618KIT_anode/R18_52189-v01.epos" "$TOOLBOX_DIR/exampleFiles/R18_52189/R18_52189-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Si/Si.RNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Si/Si.RNG" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_si/Si.RNG" "$TOOLBOX_DIR/exampleFiles/Si/Si.RNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Si/Si.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Si/Si.pos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_si/Si.pos" "$TOOLBOX_DIR/exampleFiles/Si/Si.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Si/Si.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Si/Si.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_si/Si.RRNG" "$TOOLBOX_DIR/exampleFiles/Si/Si.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Si/Si.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Si/Si.epos" && \
  ln -s "../exampleFiles_backup/external_repos/7986279/usa_denton_smith_apav_si/Si.epos" "$TOOLBOX_DIR/exampleFiles/Si/Si.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Steel before D-charging R96_50397/Steel before D-charging R96_50397.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Steel before D-charging R96_50397/Steel before D-charging R96_50397.EPOS" && \
  ln -s "../exampleFiles_backup/external_repos/5838237/Steel before D-charging R96_50397.EPOS" "$TOOLBOX_DIR/exampleFiles/Steel before D-charging R96_50397/Steel before D-charging R96_50397.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/a6c565c8-a27c-4502-93e1-96ecefa315eb.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/a6c565c8-a27c-4502-93e1-96ecefa315eb.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/a6c565c8-a27c-4502-93e1-96ecefa315eb.POS" "$TOOLBOX_DIR/exampleFiles/a6c565c8-a27c-4502-93e1-96ecefa315eb/a6c565c8-a27c-4502-93e1-96ecefa315eb.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7.POS" "$TOOLBOX_DIR/exampleFiles/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7/4eef8bf6-c36e-4d1e-95e3-e38174fcb4b7.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/aa137df8-ce0e-481e-93e6-a22cb5a34882.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/aa137df8-ce0e-481e-93e6-a22cb5a34882.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/aa137df8-ce0e-481e-93e6-a22cb5a34882.POS" "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/aa137df8-ce0e-481e-93e6-a22cb5a34882.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/aa137df8-ce0e-481e-93e6-a22cb5a34882/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5086_258504/R5086_258504_UNSM_80micron_Tip2 Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5086_258504/R5086_258504_UNSM_80micron_Tip2 Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R5086_258504_UNSM_80micron_Tip2 Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R5086_258504/R5086_258504_UNSM_80micron_Tip2 Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5076_68722/R5076_68722.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5076_68722/R5076_68722.epos" && \
  ln -s "../exampleFiles_backup/external_repos/R5076_68722.epos" "$TOOLBOX_DIR/exampleFiles/R5076_68722/R5076_68722.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/a086e6ba-ac96-40ac-b620-2acc9aa82806.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/a086e6ba-ac96-40ac-b620-2acc9aa82806.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/a086e6ba-ac96-40ac-b620-2acc9aa82806.POS" "$TOOLBOX_DIR/exampleFiles/a086e6ba-ac96-40ac-b620-2acc9aa82806/a086e6ba-ac96-40ac-b620-2acc9aa82806.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R17_94469/R17_94469 WRe Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R17_94469/R17_94469 WRe Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R17_94469 WRe Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R17_94469/R17_94469 WRe Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5100_228062/R5100_228062 W.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5100_228062/R5100_228062 W.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R5100_228062 W.pos" "$TOOLBOX_DIR/exampleFiles/R5100_228062/R5100_228062 W.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/2653fa9c-6c43-4481-94f4-10b67a5c682b.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/2653fa9c-6c43-4481-94f4-10b67a5c682b.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/2653fa9c-6c43-4481-94f4-10b67a5c682b.POS" "$TOOLBOX_DIR/exampleFiles/2653fa9c-6c43-4481-94f4-10b67a5c682b/2653fa9c-6c43-4481-94f4-10b67a5c682b.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R76_31053/R76_31053-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R76_31053/R76_31053-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/2669454/R76_31053-v01.epos" "$TOOLBOX_DIR/exampleFiles/R76_31053/R76_31053-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R17_101977/R17_101977 SiN.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R17_101977/R17_101977 SiN.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R17_101977 SiN.pos" "$TOOLBOX_DIR/exampleFiles/R17_101977/R17_101977 SiN.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/APT_TEM_dataset_selection/APT_TEM_dataset_selection.ato" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/APT_TEM_dataset_selection/APT_TEM_dataset_selection.ato" && \
  ln -s "../exampleFiles_backup/external_repos/APT_TEM_dataset_selection.ato" "$TOOLBOX_DIR/exampleFiles/APT_TEM_dataset_selection/APT_TEM_dataset_selection.ato" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_54176/R18_54176-v01.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_54176/R18_54176-v01.pos" && \
  ln -s "../exampleFiles_backup/150320UCSB_ZrO_shadow mask/R18_54176-v01.pos" "$TOOLBOX_DIR/exampleFiles/R18_54176/R18_54176-v01.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_54176/54176.rrng" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_54176/54176.rrng" && \
  ln -s "../exampleFiles_backup/150320UCSB_ZrO_shadow mask/54176.rrng" "$TOOLBOX_DIR/exampleFiles/R18_54176/54176.rrng" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R56_138794/R56_138794 YWT Top Level ROI.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R56_138794/R56_138794 YWT Top Level ROI.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R56_138794 YWT Top Level ROI.pos" "$TOOLBOX_DIR/exampleFiles/R56_138794/R56_138794 YWT Top Level ROI.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/03297260-c180-425c-8896-b0484af71750.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/03297260-c180-425c-8896-b0484af71750.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/03297260-c180-425c-8896-b0484af71750.POS" "$TOOLBOX_DIR/exampleFiles/03297260-c180-425c-8896-b0484af71750/03297260-c180-425c-8896-b0484af71750.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R18_61451/R18_61451 matchUSYD.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R18_61451/R18_61451 matchUSYD.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R18_61451 matchUSYD.pos" "$TOOLBOX_DIR/exampleFiles/R18_61451/R18_61451 matchUSYD.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/451a48eb-0d2c-4b6c-8587-ebfc832e38ff.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/451a48eb-0d2c-4b6c-8587-ebfc832e38ff.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/451a48eb-0d2c-4b6c-8587-ebfc832e38ff.POS" "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/451a48eb-0d2c-4b6c-8587-ebfc832e38ff.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/451a48eb-0d2c-4b6c-8587-ebfc832e38ff/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R6001_233217/R6001_233217CuBe.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R6001_233217/R6001_233217CuBe.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R6001_233217CuBe.pos" "$TOOLBOX_DIR/exampleFiles/R6001_233217/R6001_233217CuBe.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/4c1e71dc-9fdb-44ab-b7fb-105235775029.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/4c1e71dc-9fdb-44ab-b7fb-105235775029.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/4c1e71dc-9fdb-44ab-b7fb-105235775029.POS" "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/4c1e71dc-9fdb-44ab-b7fb-105235775029.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/4c1e71dc-9fdb-44ab-b7fb-105235775029/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/e7923146-1755-47ec-876b-ded31b78d2b6.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/e7923146-1755-47ec-876b-ded31b78d2b6.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/e7923146-1755-47ec-876b-ded31b78d2b6.POS" "$TOOLBOX_DIR/exampleFiles/e7923146-1755-47ec-876b-ded31b78d2b6/e7923146-1755-47ec-876b-ded31b78d2b6.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R5076_69132/R5076_69132-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R5076_69132/R5076_69132-v01.epos" && \
  ln -s "../exampleFiles_backup/external_repos/10677562/R5076_69132-v01.epos" "$TOOLBOX_DIR/exampleFiles/R5076_69132/R5076_69132-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_13453/R02_13453-v01_Cuboidal.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_13453/R02_13453-v01_Cuboidal.pos" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/R02_13453-v01_Cuboidal.pos" "$TOOLBOX_DIR/exampleFiles/R02_13453/R02_13453-v01_Cuboidal.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R02_13453/Rene N5.RNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R02_13453/Rene N5.RNG" && \
  ln -s "../exampleFiles_backup/Files_for_Peter_Felfer_interfacial_excess_calculation/Rene N5.RNG" "$TOOLBOX_DIR/exampleFiles/R02_13453/Rene N5.RNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19415/R04_19415-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19415/R04_19415-v01.epos" && \
  ln -s "../exampleFiles_backup/unranged/Ni/R04_19415-v01.epos" "$TOOLBOX_DIR/exampleFiles/R04_19415/R04_19415-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/98e28a97-ad0c-456a-a561-70e0a961eaf5.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/98e28a97-ad0c-456a-a561-70e0a961eaf5.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/98e28a97-ad0c-456a-a561-70e0a961eaf5.POS" "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/98e28a97-ad0c-456a-a561-70e0a961eaf5.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/98e28a97-ad0c-456a-a561-70e0a961eaf5/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/2min FeO R76_50312/2min FeO R76_50312.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/2min FeO R76_50312/2min FeO R76_50312.EPOS" && \
  ln -s "../exampleFiles_backup/external_repos/5838237/2min FeO R76_50312.EPOS" "$TOOLBOX_DIR/exampleFiles/2min FeO R76_50312/2min FeO R76_50312.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R17_98556/R17_98556.pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R17_98556/R17_98556.pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R17_98556.pos" "$TOOLBOX_DIR/exampleFiles/R17_98556/R17_98556.pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R04_19413/R04_19413-v01.epos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R04_19413/R04_19413-v01.epos" && \
  ln -s "../exampleFiles_backup/unranged/Ni/R04_19413-v01.epos" "$TOOLBOX_DIR/exampleFiles/R04_19413/R04_19413-v01.epos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/Steel after D-charging R96_50473/Steel after D-charging R96_50473.EPOS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/Steel after D-charging R96_50473/Steel after D-charging R96_50473.EPOS" && \
  ln -s "../exampleFiles_backup/external_repos/5838237/Steel after D-charging R96_50473.EPOS" "$TOOLBOX_DIR/exampleFiles/Steel after D-charging R96_50473/Steel after D-charging R96_50473.EPOS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee.POS" "$TOOLBOX_DIR/exampleFiles/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee/e89fa9d3-eef5-49e0-b0b6-9da0ca826eee.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/0d2cc974-9d9a-44ab-a779-09740f6ecb82.POS" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/0d2cc974-9d9a-44ab-a779-09740f6ecb82.POS" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/0d2cc974-9d9a-44ab-a779-09740f6ecb82.POS" "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/0d2cc974-9d9a-44ab-a779-09740f6ecb82.POS" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/range file.RRNG" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/range file.RRNG" && \
  ln -s "../exampleFiles_backup/external_repos/5534859/range file.RRNG" "$TOOLBOX_DIR/exampleFiles/0d2cc974-9d9a-44ab-a779-09740f6ecb82/range file.RRNG" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi
if [ -L "$TOOLBOX_DIR/exampleFiles/R6006_254275/R6006_254275 Cu LP .pos" ]; then
  rm "$TOOLBOX_DIR/exampleFiles/R6006_254275/R6006_254275 Cu LP .pos" && \
  ln -s "../exampleFiles_backup/external_repos/globus/R6006_254275 Cu LP .pos" "$TOOLBOX_DIR/exampleFiles/R6006_254275/R6006_254275 Cu LP .pos" && \
  fixed=$((fixed + 1)) || failed=$((failed + 1))
fi

echo "Fixed: $fixed, Failed: $failed"