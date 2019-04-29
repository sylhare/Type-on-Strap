# coding: utf-8

Gem::Specification.new do |spec|
  spec.name             = "type-on-strap"
  spec.version          = "1.2.1"
  spec.license          = ['MIT']
  spec.licenses         = ['MIT', 'OFL-1.1']
  spec.authors          = ["Sylhare","Rohan Chandra"]
  spec.email            = ["sylhare@outlook.com", "hellorohan@outlook.com"]

  spec.homepage      = "https://github.com/sylhare/Type-on-Strap"
  spec.summary          =  "A fast, simple and responsive jekyll theme template"
  spec.description      = %q{A custom Type Theme template (a free and open-source Jekyll theme). Great for blogs, easy to customize and responsive.}
  spec.extra_rdoc_files = ['README.md', 'LICENSE']

  spec.post_install_message = "Thanks for using type-on-strap!"
  spec.metadata = {
    "bug_tracker_uri"   => "https://github.com/sylhare/Type-on-Strap/issues",
    "changelog_uri"     => "https://github.com/sylhare/Type-on-Strap/releases",
    "homepage_uri"      => "https://github.com/sylhare/Type-on-Strap",
    "source_code_uri"   => "https://github.com/sylhare/Type-on-Strap",
    "wiki_uri"          => "https://github.com/sylhare/Type-on-Strap/wiki"
  }

  spec.files            = Dir.glob("**/{*,.*").select do |f|
    f.match(%r{^(assets/(js|css|fonts|data)|_(includes|layouts|sass)/|(LICENSE|README.md))})
  end

  spec.required_ruby_version = '~> 2.3'
    
  spec.add_runtime_dependency "jekyll", "~> 3.8", ">= 3.8.5"
  spec.add_runtime_dependency "jekyll-paginate", "~> 1.1"
  spec.add_runtime_dependency "jekyll-seo-tag", "~>2.6"

  spec.add_development_dependency "bundler", "~> 2.0", ">= 2.0.1"
  spec.add_development_dependency "rake", "~> 10.0"

end
