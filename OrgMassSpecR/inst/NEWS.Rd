\name{NEWS}
\title{OrgMassSpecR News}

\section{Changes in version 0.4-4 (2014-05-21)}{

    \itemize{

      \item Fixed error in ConvertPeptide. When output = "elements", the resulting element list did not include the additon of water. This caused the mass calculated by, for example, MonoisotopicMass(ConvertPeptide("SEQENCE")) to be incorrect.

    }

}

\section{Changes in version 0.4-2 (2014-03-14)}{

    \itemize{

      \item Vingette changed to R Markdown format and complied using knitr.

    }

}

\section{Changes in version 0.4 (2013-10-04)}{

    \itemize{

      \item Added trypsin.strict option to Digest function (thanks to S. Michalik).

    }

}

\section{Changes in version 0.3-18 (2013-09-04)}{

    \itemize{

      \item Function Digest fixed. It now correctly determines peptide sequences when the number of missed cleavages is > 1 (thanks to S. Michalik).

    }

}

\section{Changes in version 0.3-11 (2012-01-29)}{

    \itemize{

      \item R >= 2.14.1 required due to changes in the vignette file
      structure.

    }

}

\section{Changes in version 0.3-8 (2011-10-10)}{

  \itemize{

    \item Function ReadMspFile added.

  }

}

\section{Changes in version 0.3-6 (2011-08-26)}{

  \itemize{
  
    \item Function ConvertConcentration added.

  }

}

\section{Changes in version 0.3-2 (2011-04-08)}{

  \itemize{
  
    \item Warning message in ListFormula modified to include the
    elemental formula (helps in tracking down the source of the warning
    when called multiple times within another function).

  }

}
  

\section{Changes in version 0.3-1 (2011-03-18)}{

  \itemize{

    \item NEWS file started.

  }
  
}
