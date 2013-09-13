(when (> (length *posix-argv*) 1)
  (let ((file (second *posix-argv*)))
    (when (length file)
      (let ((geom (with-open-file (f file) (read f))))
        (format t "~a~%" (car (nth 3 geom)))))))
  
