.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    { # > figlet mascaRet
      packageStartupMessage(
"            _____                      _ 
           |  __ \\                    | |
   ___ __ _| |__) |__ _ _ __ ___   ___| |
  / __/ _` |  _  // _` | '_ ` _ \\ / _ \\ |
 | (_| (_| | | \\ \\ (_| | | | | | |  __/ |
  \\___\\__,_|_|  \\_\\__,_|_| |_| |_|\\___|_|, version ", version)

}
else
  { packageStartupMessage("Package 'caRamel' version ", version) } 
  invisible()
}
  
