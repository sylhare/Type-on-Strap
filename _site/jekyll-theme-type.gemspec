# coding: utf-8

Gem::Specification.new do |spec|
  spec.name          = "jekyll-theme-type"
  spec.version       = "1.1"
  spec.authors       = ["Rohan Chandra"]
  spec.email         = ["hellorohan@outlook.com"]

  spec.summary       = %q{A free and open-source Jekyll theme. Great for blogs and easy to customize.}
  spec.homepage      = "https://github.com/rohanchandra/type-theme"
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").select { |f| f.match(%r{^(assets|_layouts|_includes|_sass|LICENSE|README)}i) }

  spec.add_runtime_dependency "jekyll", "~> 3.4"
  spec.add_runtime_dependency "jekyll-paginate", "~> 1.1"

  spec.add_development_dependency "bundler", "~> 1.12"
  spec.add_development_dependency "rake", "~> 10.0"
end
