set(SCRIPTS_DIR Scripts)
set(SCRIPTS qi_composer.sh)
foreach(SCRIPT ${SCRIPTS})
    INSTALL( FILES ${SCRIPT} 
             DESTINATION bin 
             PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
                         OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE )
endforeach(SCRIPT)