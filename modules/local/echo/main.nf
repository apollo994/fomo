process ECHO_FILE {
    tag "${input_file.baseName}"

    input:
    path input_file

    output:
    path "*.echo.txt", emit: echoed_file

    script:
    """
    cat ${input_file} > ${input_file.baseName}.echo.txt
    """
}
