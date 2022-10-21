# coding: utf-8

Gem::Specification.new do |spec|
  spec.name          = "type-on-strap"
  spec.version       = "2.4.4"
  spec.authors       = ["Sylhare", "Rohan Chandra"]
  spec.email         = ["sylhare@outlook.com", "hellorohan@outlook.com"]

  spec.summary       = "A simple and responsive jekyll theme template"
  spec.description   = %q{A simple and responsive jekyll theme template based on type-theme. Great for blogs, easy to customize and responsive.}
  spec.homepage      = "https://github.com/sylhare/Type-on-Strap"
  spec.license       = "MIT"

  spec.rdoc_options            = ["--charset=UTF-8"]
  spec.extra_rdoc_files        = %w(README.md LICENSE)
  spec.metadata["plugin_type"] = "theme"

  spec.files                   = Dir["**/*"].select do |f|
    f.match(%r!^(assets/(js|css|fonts|data)/|_(includes|layouts|sass)/|_data/(icons_builder.yml|language.yml)|(LICENSE|README.md))!i)
  end

  spec.post_install_message =  <<~MSG
                               Thanks for using Type on strap v2+!
                               MSG

  spec.required_ruby_version   = '>= 2.7.0'

  spec.add_runtime_dependency "jekyll", ">= 3.8", "< 5.0"
  spec.add_runtime_dependency "jekyll-feed", ">= 0.15.1", "<= 0.17"
  spec.add_runtime_dependency "jekyll-paginate", "~> 1.1.0"
  spec.add_runtime_dependency "jekyll-seo-tag", ">= 2.7.1", "<= 2.8"
end
