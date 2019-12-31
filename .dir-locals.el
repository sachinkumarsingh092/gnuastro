;; This files contains Emacs Directory Local Variables.
;;
;; Emacs is an extensible, customizable, free/libre text editor.  It
;; allows specification of certain settings that will be applied to
;; all files in current directory and its subdirectories. This is
;; useful in order to automatically enforce certain coding conventions
;; for all contributors of Gnuastro, like the maximum length of lines
;; or the number of spaces to be used for indentation.
;;
;; For more information see (info "(emacs) Directory Variables")
;;
;; Copyright (C) 2015-2020 Free Software Foundation, Inc.
;;
;; Copying and distribution of this file, with or without modification,
;; are permitted in any medium without royalty provided the copyright
;; notice and this notice are preserved.  This file is offered as-is,
;; without any warranty.

((nil
  (indent-tabs-mode . nil) ;; No tabs as indentation
  (fill-column . 75))      ;; 75-character wide lines
 (c-mode
  (c-basic-offset . 2)     ;; 2 spaces of indentation
  (c-file-style . "gnu"))  ;; GNU style for braces
 (makefile-mode
  (indent-tabs-mode . t))  ;; Real TABs are important in makefiles
 )
