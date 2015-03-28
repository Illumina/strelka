
# update file with various build-time properties
# requires CONFIG_FILE SOURCE_FILE DEST_FILE

file (READ ${CONFIG_FILE} CONFIG_LINES)
STRING(REPLACE "\n" ";" CONFIG_LINES "${CONFIG_LINES}")
foreach (CONFIG_LINE ${CONFIG_LINES})
    STRING(REPLACE "\t" ";" CONFIG_LIST "${CONFIG_LINE}")
    list (GET CONFIG_LIST 0 PAIR0)
    list (GET CONFIG_LIST 1 PAIR1)
    set(${PAIR0} "${PAIR1}")
endforeach ()
configure_file(${SOURCE_FILE} ${DEST_FILE} @ONLY)
