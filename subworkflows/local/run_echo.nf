include { ECHO_FILE } from '../../modules/local/echo/main'

workflow RUN_ECHO {
    take:
    input_files

    main:
    ECHO_FILE(input_files)

    emit:
    echoed_files = ECHO_FILE.out.echoed_file
}
