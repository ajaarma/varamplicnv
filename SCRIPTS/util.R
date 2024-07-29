require("optparse")
getCommandOptionCheck = function(opt,opt_parser) {
    
    if (is.null(opt$ampCountMat)){
        print_help(opt_parser)
        stop("Amplicon count matrix R-data is missing\n",call.=FALSE)
    }
    else if (is.null(opt$bedFile)){
        print_help(opt_parser)
        stop("Amplicon bed file missing\n",call.=FALSE)
    }
    else if (is.null(opt$roi)){
        print_help(opt_parser)
        stop("ROI file is missing\n",call.=FALSE)
    }
    else if (is.null(opt$gcContent)){
        print_help(opt_parser)
        stop("GC-Content file is missing\n",call.=FALSE)
    }
    else if (is.null(opt$gender)){
        print_help(opt_parser)
        stop("Sample-Gender file is missing\n",call.=FALSE)
    }
    else if (is.null(opt$batch)){
        print_help(opt_parser)
        stop("Analysis batch path is missing\n",call.=FALSE)
    }
    else if (is.null(opt$propVarFrac)){
        print('here')
        print_help(opt_parser)
        stop("Proprtion of Variance is missing\n",call.=FALSE)
    }
    else if (is.null(opt$analysis)){
        print_help(opt_parser)
        stop("DS/AOF flag is missing\n",call.=FALSE)
    }

}


createFolders <- function(out_dir,anal_type) {
    
    out_anal_dir = paste(out_dir,"/",anal_type,"/OUTPUT/",sep="")
    out_plot_dir = paste(out_dir,"/",anal_type,"/PLOTS/",sep="")
    system(paste("mkdir -p ",out_anal_dir,sep=""))
    system(paste("mkdir -p ",out_plot_dir,sep=""))

    dir_list = list()
    dir_list = list(out_anal_dir,out_plot_dir)

    return(dir_list)

}
