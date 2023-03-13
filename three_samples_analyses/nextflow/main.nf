// Prova 

// Channel paths
files = Channel.fromPath('...')

// Read
process readFiles {

    path input_path

    input:
      file_path
    output:
      file

    script:
    """
    python read.py -i $file_path
    """
}

// Process
process processFiles {

    [ directives ]

    input:
      file
    output:
      processed_file
      
    script:
    """
    python process.py $file  
    """
}

// Write
process writeFiles {

    input:
      processed_file output_path
    output:
      < process outputs >
      
    script:
    """
    python write.py -f $processed_file -o $output_path
    """
}

// Run
workflow {
  files | readFiles | processFiles | writeFiles }
}




// process < name > {
// 
//   [ directives ]
// 
//   input:
//     < process inputs >
// 
//   output:
//     < process outputs >
// 
//   when:
//     < condition >
// 
//   [script|shell|exec]:
//     < user script to be executed >
// 
// }