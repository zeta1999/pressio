
#ifndef SOLVERS_LINE_SEARCH_TAGS_HPP_
#define SOLVERS_LINE_SEARCH_TAGS_HPP_

namespace rompp{ namespace solvers{ namespace iterative{ namespace gn{

/* if you add more here, rememeber to change also the metafunction:
 * is_non_default_line_search_tag inside meta
 */

struct noLineSearch{};
struct ArmijoLineSearch{};

}}}}//end namespace rompp::solvers::iterative::gn
#endif
