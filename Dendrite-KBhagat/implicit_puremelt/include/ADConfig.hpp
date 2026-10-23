#pragma once

namespace AD{
    // automatic differentiation setup
    constexpr Differentiation::AD::NumberTypes ADTypeCode = Differentiation::AD::NumberTypes::sacado_dfad;
    using ADHelper = Differentiation::AD::ResidualLinearization<ADTypeCode, double>;
    using ADNumber = typename ADHelper::ad_type;
}