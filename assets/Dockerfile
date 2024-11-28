FROM sylhare/jekyll:latest
LABEL maintainer="sylhare"
LABEL image="sylhare/type-on-strap"

# Create Type-on-strap Gemfile
RUN echo "source \"https://rubygems.org\"" >> Gemfile
RUN echo "gem 'type-on-strap', '>= 2.4.11', '< 3.0'" >> Gemfile
RUN echo "Adding the Gemfile"
RUN cat Gemfile

# Install the theme
RUN bundle update
RUN bundle install

# Make it accessible from outside
EXPOSE 4000

