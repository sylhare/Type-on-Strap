# coding: utf-8

Gem::Specification.new do |spec|
  spec.name          = "type-on-strap"
  spec.version       = "2.0.0"
  spec.authors       = ["Sylhare","Rohan Chandra"]
  spec.email         = ["sylhare@outlook.com", "hellorohan@outlook.com"]

  spec.summary       =  "A simple and responsive jekyll theme template"
  spec.description   = %q{A simple and responsive jekyll theme template based on type-theme. Great for blogs, easy to customize and responsive.}
  spec.homepage      = "https://github.com/sylhare/Type-on-Strap"
  spec.license       = "MIT"

  spec.rdoc_options            = ["--charset=UTF-8"]
  spec.extra_rdoc_files        = %w(README.md LICENSE)
  spec.metadata["plugin_type"] = "theme"

  spec.files                   = `git ls-files -z`.split("\x0").select do |f|
    f.match(%r!^(assets/(js|css|fonts)/|_(includes|layouts|sass)/|_data/language.yml|(LICENSE|README)((\.(txt|md|markdown)|$)))!i)
  end

    spec.post_install_message =  <<~MSG
                                    ----------------------------------------------------------

                                    Type on strap v2+ is using Jekyll 4.0:

                                      * Please make sure you have updated your _config.yml.

                                      * Use _data/ for social and language customization

                                      * Vist https://github.com/sylhare/Type-on-Strap
                                      for more info.

                                    ----------------------------------------------------------
                                  MSG

  spec.required_ruby_version   = '>= 2.4.0'
    
  spec.add_runtime_dependency "jekyll", ">= 3.5", "< 5.0"
  spec.add_runtime_dependency "jekyll-feed", "~> 0.9"
  spec.add_runtime_dependency "jekyll-paginate", "~> 1.1"
  spec.add_runtime_dependency "jekyll-seo-tag", "~>2.6"

  spec.add_development_dependency "bundler", ">= 2.0", "< 2.1.0"
end
