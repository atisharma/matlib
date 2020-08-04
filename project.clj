;; Copyright (c) 2020 A S Sharma. All rights reserved.
;;
;; The use and distribution terms for this software are covered by the
;; Eclipse Public License 2.0 or later which can be found in the file LICENSE
;; at the root of this distribution.
;; By using this software in any fashion, you are agreeing to be bound by
;; the terms of this license.
;; You must not remove this notice, or any other, from this software.

(defproject matlib "0.1.8"
  :description "A Clojure library of optimisation and control theory tools and convenience functions based on Neanderthal."
  :license {:name "EPL-2.0 OR GPL-2.0-or-later WITH Classpath-exception-2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  :scm {:name "git"
        :url "https://github.com/atisharma/matlib"}
  :url "https://agalmic.ltd"
  ; jvm-opts required by Neanderthal on JDK > 8.
  :jvm-opts ^:replace ["--add-opens=java.base/jdk.internal.ref=ALL-UNNAMED"]
  :dependencies [[org.clojure/clojure "1.10.1"]
                 [uncomplicate/neanderthal "0.34.0"]
                 ; No need to specify slf4j-api, as it’s required by logback.
                 ; These loggers are here mainly to keep Neanderthal quiet.
                 [org.clojure/tools.logging "0.6.0"]
                 [ch.qos.logback/logback-classic "1.2.3"]]
  ; lein with-env-vars repl
  ; Is this working? Specified in $PATH in .bashrc.
  :env-vars {:LD_LIBRARY_PATH "/opt/intel/mkl/lib:/opt/intel/lib:/opt/intel/mkl/lib/intel64:/opt/intel/mkl/lib/intel64_lin"}
  :repositories [["local" "~/.m2/repository"
                  :creds :gpg]]
  :signing {:gpg-key "ati@agalmic.ltd"}
  :codox {:metadata {:doc/format :markdown}}
  :repl-options {:init-ns matlib.core}
  :profiles {:repl {:dependencies [[agalmic/plot "0.1.5"]]}})
