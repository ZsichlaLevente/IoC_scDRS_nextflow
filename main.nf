#!/usr/bin/env nextflow

process sayHello {
    container 'marcoschinas/scdrs:1.0.2'

    input:
    val name

    output:
    stdout

    script:
    """
    echo "Hello $name from scDRS!"
    """
}

workflow {
    sayHello('Levente')
}
