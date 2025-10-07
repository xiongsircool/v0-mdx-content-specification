import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { CalendarDays, Clock, User, ExternalLink, AlertTriangle, MapPin, Users, Video } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import Image from "next/image"
import type { ContentItem } from "@/lib/server-markdown-loader"

interface ContentCardProps {
  item: ContentItem
  type: "announcement" | "post" | "resource" | "meeting"
}

export function ContentCard({ item, type }: ContentCardProps) {
  const getDetailUrl = () => {
    switch (type) {
      case "announcement":
        return `/announcements/${item.slug}`
      case "post":
        return `/posts/${item.slug}`
      case "resource":
        return `/learn/resources/${item.slug}`
      case "meeting":
        return `/meetings/${item.slug}`
      default:
        return "#"
    }
  }

  const getExcerptOrSummary = () => {
    return item.excerpt || item.summary || ""
  }

  const getTags = () => {
    return item.tags || []
  }

  const getDateInfo = () => {
    if (type === "meeting" && item.date) {
      return {
        icon: CalendarDays,
        text: new Date(item.date).toLocaleDateString('zh-CN'),
        label: "会议时间"
      }
    }

    if (item.publishedAt) {
      return {
        icon: CalendarDays,
        text: new Date(item.publishedAt).toLocaleDateString('zh-CN'),
        label: "发布时间"
      }
    }

    if (item.updatedAt) {
      return {
        icon: Clock,
        text: new Date(item.updatedAt).toLocaleDateString('zh-CN'),
        label: "更新时间"
      }
    }

    return null
  }

  const dateInfo = getDateInfo()

  return (
    <Card className="transition-all hover:shadow-md">
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1">
            <h3 className="text-lg font-semibold leading-tight mb-2">
              <Link href={getDetailUrl()} className="hover:text-primary transition-colors">
                {item.title}
              </Link>
            </h3>

            {getExcerptOrSummary() && (
              <p className="text-sm text-muted-foreground line-clamp-2 mb-3">
                {getExcerptOrSummary()}
              </p>
            )}

            <div className="flex flex-wrap gap-2">
              {getTags().map(tag => (
                <Badge key={tag} variant="secondary" className="text-xs">
                  {tag}
                </Badge>
              ))}
            </div>
          </div>

          {item.coverImage && (
            <div className="w-20 h-20 flex-shrink-0">
              <Image
                src={item.coverImage}
                alt={item.title}
                width={80}
                height={80}
                className="w-full h-full object-cover rounded-md"
              />
            </div>
          )}
        </div>

        {dateInfo && (
          <div className="flex items-center gap-2 text-sm text-muted-foreground mt-3">
            <dateInfo.icon className="h-4 w-4" />
            <span>{dateInfo.text}</span>
          </div>
        )}

        {/* Meeting specific info */}
        {type === "meeting" && (
          <div className="flex flex-wrap gap-3 text-sm text-muted-foreground mt-2">
            {item.location && (
              <div className="flex items-center gap-1">
                <MapPin className="h-4 w-4" />
                <span>{item.location}</span>
              </div>
            )}
            {item.eventType && (
              <Badge variant="outline" className="text-xs">
                {item.eventType}
              </Badge>
            )}
            {item.status && (
              <Badge
                variant={item.status === "即将开始" ? "default" : "secondary"}
                className="text-xs"
              >
                {item.status}
              </Badge>
            )}
          </div>
        )}

        {/* Resource specific info */}
        {type === "resource" && (
          <div className="flex flex-wrap gap-3 text-sm text-muted-foreground mt-2">
            {item.type && (
              <Badge variant="outline" className="text-xs">
                {item.type}
              </Badge>
            )}
            {item.level && (
              <Badge variant="secondary" className="text-xs">
                {item.level}
              </Badge>
            )}
            {item.duration && (
              <span>{item.duration}</span>
            )}
          </div>
        )}

        {/* Post specific info */}
        {type === "post" && item.authors && item.authors.length > 0 && (
          <div className="flex items-center gap-2 text-sm text-muted-foreground mt-2">
            <User className="h-4 w-4" />
            <span>{item.authors.join(", ")}</span>
          </div>
        )}

        {/* Announcement expiration warning */}
        {type === "announcement" && item.isExpired && (
          <div className="flex items-center gap-2 text-amber-600 text-sm mt-2">
            <AlertTriangle className="h-4 w-4" />
            <span>此公告已过期</span>
          </div>
        )}
      </CardHeader>

      <CardContent className="pt-0">
        <div className="flex items-center justify-between">
          <Button variant="ghost" size="sm" asChild>
            <Link href={getDetailUrl()}>
              查看详情
              <ExternalLink className="ml-2 h-4 w-4" />
            </Link>
          </Button>

          {/* Meeting registration link */}
          {type === "meeting" && item.meetingUrl && (
            <Button variant="outline" size="sm" asChild>
              <Link href={item.meetingUrl} target="_blank" rel="noopener noreferrer">
                <Video className="mr-2 h-4 w-4" />
                加入会议
              </Link>
            </Button>
          )}

          {/* Resource external link */}
          {type === "resource" && item.externalUrl && (
            <Button variant="outline" size="sm" asChild>
              <Link href={item.externalUrl} target="_blank" rel="noopener noreferrer">
                <ExternalLink className="mr-2 h-4 w-4" />
                访问资源
              </Link>
            </Button>
          )}
        </div>
      </CardContent>
    </Card>
  )
}